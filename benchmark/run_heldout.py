"""Run the frozen cavity-p95 method against the preregistered held-out cohort.

Order matters and is enforced here: the two pre-specified filters are applied
and reported *before* any AUROC is computed, so the eligible set cannot be
chosen after seeing which entries score well.

Labels come from a different crystal form than the trajectory, so they are
carried across by UniProt residue number: SIFTS gives the author->UniProt offset
for the labelled chain, and ATLAS's own UnP_num column gives UniProt->residue
index for the trajectory. No sequence alignment is involved.

See PREREGISTRATION.md. Prediction: mean AUROC >= 0.70.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time
import zipfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))

import mdtraj as md
from build_cohort import labels_to_uniprot, sifts_map, uniprot_to_index
from cryptic_dispersion.cavity import residue_cavity_volume
from cryptobench_atlas import auroc

STRIDE = 5
MIN_LABELS = 5
JOIN_Z = -2.0          # labels must be more spatially clustered than chance


def spatial_z(top, mask, rng, n=400):
    ca = top.topology.select("name CA")
    xyz = top.xyz[0, ca, :] * 10.0
    idx = np.flatnonzero(mask[:len(ca)])
    if len(idx) < 3 or len(idx) >= len(ca):
        return float("nan")
    def spread(ii):
        p = xyz[ii]
        return float(np.mean(np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)))
    obs = spread(idx)
    null = [spread(rng.choice(len(ca), len(idx), replace=False)) for _ in range(n)]
    return (obs - np.mean(null)) / np.std(null)


def main() -> int:
    cohort = json.loads((HERE / "heldout.json").read_text())
    labels = json.loads((HERE / "cb_labels.json").read_text())
    sifts = sifts_map(HERE / "sifts.tsv")
    zips = HERE / "heldout_cohort"
    work = zips / "_extracted"
    out_path = HERE / "heldout_result.json"
    rng = np.random.default_rng(0)

    done = {}
    if out_path.exists():
        done = {d["atlas_entry"]: d for d in json.loads(out_path.read_text())}
        print(f"resuming: {len(done)} entries scored", flush=True)

    # ---- filters, applied and reported before any score ----
    print("PRE-SPECIFIED FILTERS (applied before scoring)\n", flush=True)
    print(f"{'entry':<9} {'labels':>6} {'mapped':>7} {'join z':>7}  status", flush=True)
    eligible = []
    for e in cohort:
        tag = e["atlas_entry"]
        z = zips / f"{tag}.zip"
        if not z.exists():
            print(f"{tag:<9} {'':>6} {'':>7} {'':>7}  SKIP archive not downloaded", flush=True)
            continue
        d = work / tag
        if not (d / f"{tag}.pdb").exists():
            d.mkdir(parents=True, exist_ok=True)
            zipfile.ZipFile(z).extractall(d)
        if e["apo_guard"] == "known_holo":
            print(f"{tag:<9} {'':>6} {'':>7} {'':>7}  EXCLUDED known holo", flush=True)
            continue
        top = md.load(str(d / f"{tag}.pdb"))
        lpdb, lch = e["label_entry"].split("_")
        up_nums = labels_to_uniprot(labels[lpdb.lower()], lch, lpdb, e["uniprot"], sifts)
        u2i = uniprot_to_index(d / f"{tag}_corresp.tsv")
        idx = sorted({u2i[u] for u in up_nums if u in u2i and u2i[u] < top.n_residues})
        mask = np.zeros(top.n_residues, bool); mask[idx] = True
        zz = spatial_z(top, mask, rng)
        if len(idx) < MIN_LABELS:
            st = f"EXCLUDED only {len(idx)} labels mapped"
        elif not np.isfinite(zz) or zz > JOIN_Z:
            st = f"EXCLUDED join not clustered"
        else:
            st = "eligible"
            eligible.append((tag, d, mask, e))
        print(f"{tag:<9} {len(up_nums):>6} {len(idx):>7} {zz:>7.1f}  {st}", flush=True)

    print(f"\neligible after filters: {len(eligible)}/{len(cohort)}"
          f"  ({len({e['uniprot'] for _,_,_,e in eligible})} UniProt clusters)\n", flush=True)

    # ---- frozen scoring ----
    print("FROZEN METHOD: cavity p95, probe 1.4, enclosure 3, span 8, 3 replicas\n", flush=True)
    print(f"{'entry':<9} {'res':>5} {'lab':>4} {'CAVITY':>7} {'sasa':>6} {'rmsf':>6} {'sec':>5}", flush=True)
    print("-" * 50, flush=True)
    for tag, d, mask, e in eligible:
        if tag in done:
            continue
        t0 = time.time()
        try:
            reps = [md.load(str(d / f"{tag}_R{i}.xtc"), top=str(d / f"{tag}.pdb"))[::STRIDE]
                    for i in (1, 2, 3)]
            t = reps[0]
            heavy = [a.index for a in t.topology.atoms if a.element.symbol != "H"]
            el = [t.topology.atom(i).element.symbol for i in heavy]
            res = np.array([t.topology.atom(i).residue.index for i in heavy])
            ser = np.concatenate([
                np.array([residue_cavity_volume(r.xyz[f, heavy, :] * 10.0, el, res,
                                                t.n_residues) for f in range(r.n_frames)])
                for r in reps])
            cav = np.percentile(ser, 95, axis=0)
            sasa = md.shrake_rupley(t, mode="residue", n_sphere_points=96).mean(0)
            rmsf = md.rmsf(t, t, 0, atom_indices=t.topology.select("name CA"))
            def a(v):
                o = np.argsort(-v)
                return auroc(v[o], mask[:len(v)][o])
            rec = dict(atlas_entry=tag, label_entry=e["label_entry"],
                       uniprot=e["uniprot"], n_residues=int(t.n_residues),
                       n_labels=int(mask.sum()), cavity_auroc=a(cav),
                       sasa_auroc=a(sasa), rmsf_auroc=a(rmsf),
                       seconds=time.time() - t0)
        except Exception as exc:                                   # noqa: BLE001
            print(f"{tag:<9} FAILED {type(exc).__name__}: {exc}", flush=True)
            continue
        done[tag] = rec
        out_path.write_text(json.dumps(list(done.values()), indent=2))
        print(f"{tag:<9} {rec['n_residues']:>5} {rec['n_labels']:>4} "
              f"{rec['cavity_auroc']:>7.3f} {rec['sasa_auroc']:>6.3f} "
              f"{rec['rmsf_auroc']:>6.3f} {rec['seconds']:>5.0f}", flush=True)
    print("\nSCORING COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
