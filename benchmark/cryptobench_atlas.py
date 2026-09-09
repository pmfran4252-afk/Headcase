"""Score real ATLAS molecular dynamics against CryptoBench cryptic-site labels.

This is the first time the pipeline meets a protein.  Everything before it ran
on generators, and a generator cannot tell you whether a channel is measuring
the thing it was designed to measure or an artefact of how the data was made.

Two design decisions carry most of the weight.

**The three channels read three different measurements.**  Burial (a CA contact
count), openness (Shrake-Rupley solvent-accessible surface area) and amide
environment (backbone N contacts plus hydrogen bonds) are separate geometric
quantities computed independently from the coordinates.  They are correlated
through the real conformational state, which is correct and is the premise of
consensus scoring -- but none is an algebraic transform of another, so their
agreement is evidence rather than arithmetic.  That is the distinction the
synthetic demo originally failed.

**Replicas are scored separately and averaged, not concatenated.**  ATLAS ships
three independent 100 ns runs.  Concatenating them into one 300 ns series, as
the project README suggests, triples the frame count but inserts two
discontinuities that the HMM reads as transitions and that inflate the
integrated autocorrelation time.  Pooling helps counting statistics and corrupts
every time-series estimate, so each replica is scored on its own and the
per-residue channel values are averaged across the three.

Results are written after every entry rather than at the end, and an existing
result file is read back on startup so a re-run resumes.  A previous version
buffered everything and wrote once at the end; it was killed partway and the
entire cohort was lost.

Labels: CryptoBench (Skrhak et al., Bioinformatics 2025).
Trajectories: ATLAS (Vander Meersche et al., NAR 2024), CC-BY-NC.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
import time
import zipfile
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from cryptic_dispersion.exchange import exchange_contrast, fit_two_state
from cryptic_dispersion.hdx import AmideEnvironment, breathing_anomaly
from cryptic_dispersion.observability import (
    ABSENT, BLIND, ChannelStatus, REFERENCE_K_EX, REFERENCE_POPULATION,
    ResidueStatus, breathing_status, dispersion_floor, dispersion_status,
    residue_outcome, summarise, tail_status,
)
from cryptic_dispersion.score import rank_sites
from cryptic_dispersion.tails import latent_openness

SASA_SPHERE_POINTS = 96      # r = 0.9988 against the 960-point reference, 4.3x faster
CONTACT_CUTOFF_NM = 1.00     # CA-CA burial cutoff
EXCLUDE_WINDOW = 2           # chain neighbours contribute burial that is not tertiary
DT_PS = 100.0                # ATLAS saves every 100 ps


@dataclass
class EntryResult:
    tag: str
    n_residues: int
    n_labelled: int
    label_match_rate: float
    top_k_precision: float
    base_rate: float
    enrichment: float
    auroc: float
    absent: int
    seconds: float = 0.0
    blind: dict = field(default_factory=dict)
    note: str = ""

    def as_dict(self) -> dict:
        return dict(self.__dict__)


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

def load_corresp(path: pathlib.Path) -> dict:
    """ATLAS renumbers residues from 1.  Map that back to author numbering.

    Returns ``{author_resnum: residue_index}`` for both the PDB and CIF columns,
    because CryptoBench's numbering convention is not documented in the label
    file and has to be established by which one matches.
    """
    rows = [l.split("\t") for l in path.read_text().splitlines() if l.strip()]
    head, rows = rows[0], rows[1:]
    ipdb, icif = head.index("PDB_num"), head.index("CIF_num")
    by_pdb, by_cif = {}, {}
    for idx, r in enumerate(rows):
        for col, target in ((ipdb, by_pdb), (icif, by_cif)):
            if len(r) > col:
                try:
                    target[int(r[col])] = idx
                except ValueError:
                    pass
    return {"pdb": by_pdb, "cif": by_cif}


def load_entry(zip_path: pathlib.Path, work: pathlib.Path):
    import mdtraj as md
    tag = zip_path.stem
    out = work / tag
    if not (out / f"{tag}.pdb").exists():
        out.mkdir(parents=True, exist_ok=True)
        zipfile.ZipFile(zip_path).extractall(out)
    top = str(out / f"{tag}.pdb")
    reps = [md.load(str(out / f"{tag}_R{i}.xtc"), top=top) for i in (1, 2, 3)]
    return tag, reps, load_corresp(out / f"{tag}_corresp.tsv")


# --------------------------------------------------------------------------
# Observables: three separate geometric measurements
# --------------------------------------------------------------------------

def burial(traj) -> np.ndarray:
    """Per-residue CA contact number, negated so it grows as the site opens."""
    ca = traj.topology.select("name CA")
    xyz = traj.xyz[:, ca, :].astype(np.float64)
    n = xyz.shape[1]
    seq = np.arange(n)
    near = np.abs(seq[:, None] - seq[None, :]) <= EXCLUDE_WINDOW
    out = np.empty((traj.n_frames, n))
    for s in range(0, traj.n_frames, 128):
        e = min(s + 128, traj.n_frames)
        d = np.linalg.norm(xyz[s:e, :, None, :] - xyz[s:e, None, :, :], axis=-1)
        out[s:e] = ((d < CONTACT_CUTOFF_NM) & ~near).sum(-1)
    return -out


def exposure(traj) -> np.ndarray:
    """Per-residue SASA -- a geometric quantity independent of the contact count,
    and the natural openness proxy when no pocket-detection head is available."""
    import mdtraj as md
    return md.shrake_rupley(traj, mode="residue",
                            n_sphere_points=SASA_SPHERE_POINTS).astype(np.float64)


def amide_env(traj):
    """Backbone amide burial and hydrogen bonding, per frame, in angstroms."""
    top = traj.topology
    n_idx, h_idx, res_of = [], [], []
    for r in top.residues:
        n = next((a.index for a in r.atoms if a.name == "N"), None)
        h = next((a.index for a in r.atoms if a.name in ("H", "HN", "H1")), None)
        if n is not None and h is not None:
            n_idx.append(n); h_idx.append(h); res_of.append(r.index)
    heavy = np.array([a.index for a in top.atoms if a.element.symbol != "H"])
    heavy_res = np.array([top.atom(i).residue.index for i in heavy])
    acc_pos = np.flatnonzero(
        np.array([top.atom(i).element.symbol in ("O", "N") for i in heavy]))

    xyz = traj.xyz * 10.0
    hx, nx, ax = xyz[:, heavy, :], xyz[:, n_idx, :], xyz[:, h_idx, :]
    t, r = traj.n_frames, len(n_idx)
    contacts, hbonds = np.zeros((t, r)), np.zeros((t, r))
    for i in range(r):
        far = np.abs(heavy_res - res_of[i]) > EXCLUDE_WINDOW
        sel = acc_pos[np.abs(heavy_res[acc_pos] - res_of[i]) > EXCLUDE_WINDOW]
        for s in range(0, t, 128):
            e = min(s + 128, t)
            d2 = ((hx[s:e][:, far, :] - nx[s:e, i, None, :]) ** 2).sum(-1)
            contacts[s:e, i] = (d2 < 6.5 ** 2).sum(1)
            dh2 = ((hx[s:e][:, sel, :] - ax[s:e, i, None, :]) ** 2).sum(-1)
            hbonds[s:e, i] = (dh2 < 2.4 ** 2).sum(1)
    return AmideEnvironment(contacts, np.minimum(hbonds, 1.0)), np.array(res_of)


# --------------------------------------------------------------------------
# Scoring one replica
# --------------------------------------------------------------------------

def score_replica(traj, dt_ps: float = DT_PS):
    """Score one replica, skipping any channel the trajectory cannot support.

    The dispersion channel is gated on its floor *before* the fit rather than
    after.  At ATLAS geometry -- 100 ns at 100 ps framing -- a 5% excited state
    needs k_ex >= 8.4e8 s^-1 to be fitted, three orders above the 1e6 s^-1 at
    which cryptic sites open, so the channel is blind for every residue of every
    protein in this cohort.  Fitting 276-500 Baum-Welch HMMs per replica to
    populate a channel that provably cannot speak was the dominant cost of this
    harness and bought nothing; the floor is arithmetic and can be checked first.
    """
    obs, openv = burial(traj), exposure(traj)
    env, amide_res = amide_env(traj)
    n = obs.shape[1]
    total_ns = traj.n_frames * dt_ps / 1000.0

    breath = np.full(n, np.nan)
    breath[amide_res] = breathing_anomaly(env)

    d_floor = dispersion_floor(total_ns, REFERENCE_POPULATION)
    blind = d_floor > REFERENCE_K_EX
    blind_status = ChannelStatus(
        "dispersion", BLIND, d_floor, "k_ex s^-1",
        f"a {REFERENCE_POPULATION:.0%} state needs k_ex >= {d_floor:.1e} s^-1 in "
        f"{total_ns:.0f} ns; cryptic opening runs near {REFERENCE_K_EX:.0e} s^-1")

    disp, tail = np.full(n, np.nan), np.full(n, np.nan)
    fits, statuses = [], []
    for i in range(n):
        if blind:
            fits.append(None)
            ds = blind_status
        else:
            f = fit_two_state(obs[:, i], dt_ps=dt_ps)
            fits.append(f)
            if f.reliable:
                disp[i] = exchange_contrast(f)
            ds = dispersion_status(f, total_ns)
        tail[i] = latent_openness(openv[:, i], target_population=1e-3)
        ch = {"dispersion": ds,
              "tail": tail_status(openv[:, i], total_ns),
              "breathing": breathing_status(
                  env.n_contacts[:, min(i, env.n_contacts.shape[1] - 1)])}
        statuses.append(ResidueStatus(i, residue_outcome(ch), ch))
    return disp, tail, breath, fits, statuses


# --------------------------------------------------------------------------
# Labels and metrics
# --------------------------------------------------------------------------

def map_labels(labels: list, chain: str, corresp: dict, n_res: int):
    """Project CryptoBench labels onto trajectory residue indices.

    CryptoBench writes ``<chain>_<resnum>`` but does not document whether the
    number is author or CIF numbering, so both are tried and the one that maps
    more labels wins.  The match rate is returned and reported: a label set that
    lands on the wrong residues would produce a confident, meaningless
    enrichment, and silence about the join is how that happens.
    """
    wanted = [int(x.split("_")[1]) for x in labels
              if x.split("_")[0] == chain and x.split("_")[1].lstrip("-").isdigit()]
    if not wanted:
        return np.zeros(n_res, bool), 0.0, "no labels on this chain"
    best, best_kind, best_hits = [], "", -1
    for kind in ("pdb", "cif"):
        m = corresp[kind]
        idx = [m[w] for w in wanted if w in m and m[w] < n_res]
        if len(idx) > best_hits:
            best, best_kind, best_hits = idx, kind, len(idx)
    mask = np.zeros(n_res, bool)
    mask[best] = True
    return mask, best_hits / len(wanted), f"{best_kind} numbering"


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos, neg = scores[labels], scores[~labels]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    order = np.argsort(np.concatenate([pos, neg]))
    ranks = np.empty(order.size)
    ranks[order] = np.arange(1, order.size + 1)
    return float((ranks[:pos.size].sum() - pos.size * (pos.size + 1) / 2)
                 / (pos.size * neg.size))


def run_entry(zip_path, work, labels_by_pdb, n_replicas: int = 3) -> Optional[EntryResult]:
    t0 = time.time()
    tag, reps, corresp = load_entry(zip_path, work)
    pdb_id, chain = tag.split("_")
    labels = labels_by_pdb.get(pdb_id.lower())
    if not labels:
        return None
    reps = reps[:n_replicas]
    n = reps[0].n_residues

    per = []
    for j, r in enumerate(reps, 1):
        per.append(score_replica(r))
        print(f"    {tag} replica {j}/{len(reps)} done ({time.time()-t0:.0f}s)", flush=True)
    with np.errstate(invalid="ignore"):
        disp = np.nanmean([p[0] for p in per], axis=0)
        tail = np.nanmean([p[1] for p in per], axis=0)
        breath = np.nanmean([p[2] for p in per], axis=0)
    ranked = rank_sites(disp, tail, breath, fits=per[0][3], statuses=per[0][4])

    mask, rate, note = map_labels(labels, chain, corresp, n)
    summary = summarise(per[0][4])
    if mask.sum() == 0:
        return EntryResult(tag, n, 0, rate, float("nan"), float("nan"),
                           float("nan"), float("nan"),
                           summary["counts"][ABSENT], time.time() - t0,
                           summary["blind_by_channel"], note)
    order = [s.residue for s in ranked]
    k = int(mask.sum())
    hits = sum(1 for r in order[:k] if mask[r])
    base = k / n
    return EntryResult(
        tag=tag, n_residues=n, n_labelled=k, label_match_rate=rate,
        top_k_precision=hits / k, base_rate=base,
        enrichment=(hits / k) / base if base else float("nan"),
        auroc=auroc(np.array([s.score for s in ranked]),
                    np.array([mask[r] for r in order])),
        absent=summary["counts"][ABSENT], seconds=time.time() - t0,
        blind=summary["blind_by_channel"], note=note)


def main(argv: list) -> int:
    here = pathlib.Path(__file__).resolve().parent
    cohort = pathlib.Path(argv[1]) if len(argv) > 1 else here / "atlas_cohort"
    labels_path = pathlib.Path(argv[2]) if len(argv) > 2 else here / "cb_labels.json"
    work = cohort / "_extracted"
    n_replicas = int(os.environ.get("HEADCASE_REPLICAS", "3"))
    out_path = here / f"cryptobench_atlas_result_r{n_replicas}.json"
    labels_by_pdb = json.loads(labels_path.read_text())

    done = {}
    if out_path.exists():
        for d in json.loads(out_path.read_text()):
            done[d["tag"]] = d
        print(f"resuming: {len(done)} entries already scored", flush=True)

    zips = sorted(cohort.glob("*.zip"))
    print(f"CryptoBench x ATLAS: {len(zips)} entries, {n_replicas} replica(s) each\n", flush=True)
    print(f"{'entry':<9} {'res':>4} {'lab':>4} {'join':>5} {'top-k':>6} {'base':>6} "
          f"{'enrich':>7} {'auroc':>6} {'sec':>5}  blind channels", flush=True)
    print("-" * 88, flush=True)

    for z in zips:
        if z.stem in done:
            d = done[z.stem]
            print(f"{d['tag']:<9} (cached) enrich {d['enrichment']:.2f} "
                  f"auroc {d['auroc']:.2f}", flush=True)
            continue
        try:
            r = run_entry(z, work, labels_by_pdb, n_replicas)
        except Exception as exc:                                  # noqa: BLE001
            print(f"{z.stem:<9} FAILED {type(exc).__name__}: {exc}", flush=True)
            continue
        if r is None:
            print(f"{z.stem:<9} no labels for this pdb id", flush=True)
            continue
        done[r.tag] = r.as_dict()
        # Write after every entry. A previous version wrote once at the end and
        # lost the whole cohort when it was killed partway.
        out_path.write_text(json.dumps(list(done.values()), indent=2))
        bl = ", ".join(f"{k[:4]} {v}/{r.n_residues}" for k, v in sorted(r.blind.items()))
        print(f"{r.tag:<9} {r.n_residues:>4} {r.n_labelled:>4} {r.label_match_rate:>5.0%} "
              f"{r.top_k_precision:>6.2f} {r.base_rate:>6.2f} {r.enrichment:>7.2f} "
              f"{r.auroc:>6.2f} {r.seconds:>5.0f}  {bl}", flush=True)

    ok = [d for d in done.values() if d.get("auroc") == d.get("auroc")
          and d.get("auroc") is not None]
    if ok:
        e = [d["enrichment"] for d in ok]
        a = [d["auroc"] for d in ok]
        print(f"\ncohort: {len(ok)} scored entries", flush=True)
        print(f"  median enrichment over base rate : {np.median(e):.2f}x", flush=True)
        print(f"  median AUROC                     : {np.median(a):.3f}", flush=True)
        print(f"  entries with AUROC > 0.5         : {sum(x > 0.5 for x in a)}/{len(a)}", flush=True)
        print(f"  median label join rate           : "
              f"{np.median([d['label_match_rate'] for d in ok]):.0%}", flush=True)
        print(f"\nwrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
