"""Before/after: the same cavity detector on unbiased vs cosolvent trajectories.

The comparison is only meaningful if nothing changes except the input, so this
imports the same `residue_cavity_volume`, the same p95 aggregation, and the same
label mapping the benchmark used. No parameter is re-tuned for cosolvent data.

Topology handling is deliberately explicit. The cosolvent trajectory contains
protein atoms only, and the PDB written beside it must describe exactly those
atoms. A mismatch here is silent until analysis: the first run wrote the full
13,391-atom solvated topology next to a 610-atom trajectory, which cannot be
opened together at all. Residue names are also truncated to three characters by
the PDB format, so `BENZ` is stored as `BEN` -- selecting on the name is fragile
and `topology.select("protein")` is used instead.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent / "benchmark"))

import mdtraj as md
from align_labels import transfer
from cryptic_dispersion.cavity import residue_cavity_volume
from cryptobench_atlas import auroc, load_corresp, map_labels

STRIDE = 2


def protein_topology(top_pdb: pathlib.Path, out: pathlib.Path) -> pathlib.Path:
    """Write a topology matching the protein-only trajectory, if not already."""
    full = md.load(str(top_pdb))
    sel = full.topology.select("protein")
    if full.n_atoms == len(sel):
        return top_pdb
    full.atom_slice(sel).save_pdb(str(out))
    return out


def cavity_p95(traj) -> np.ndarray:
    heavy = [a.index for a in traj.topology.atoms if a.element.symbol != "H"]
    el = [traj.topology.atom(i).element.symbol for i in heavy]
    res = np.array([traj.topology.atom(i).residue.index for i in heavy])
    ser = np.array([residue_cavity_volume(traj.xyz[f, heavy, :] * 10.0, el, res,
                                          traj.n_residues)
                    for f in range(traj.n_frames)])
    return np.percentile(ser, 95, axis=0), ser


def compare(tag, extracted, cosolv_dir, labels, label_entry):
    d = pathlib.Path(extracted) / tag
    cos = pathlib.Path(cosolv_dir)
    lpdb, lch = label_entry.split("_")

    # labels, mapped exactly as the benchmark maps them
    ref = md.load(str(d / f"{tag}.pdb"))
    lab = labels[lpdb.lower()]
    corr = load_corresp(d / f"{tag}_corresp.tsv")
    mask, rate, note = map_labels(lab, lch, corr, ref.n_residues)
    if mask.sum() < 3:
        want = {int(x.split("_")[1]) for x in lab
                if x.split("_")[0] == lch and x.split("_")[1].lstrip("-").isdigit()}
        hits, ident, _ = transfer(lpdb, lch, d / f"{tag}_corresp.tsv", want,
                                  HERE.parent / "benchmark" / "pdbcache")
        mask = np.zeros(ref.n_residues, bool)
        for h in hits:
            if h < ref.n_residues:
                mask[h] = True

    out = {"tag": tag, "n_labels": int(mask.sum())}
    # BEFORE: unbiased ATLAS
    before = md.load(str(d / f"{tag}_R1.xtc"), top=str(d / f"{tag}.pdb"))[::10]
    cb, _ = cavity_p95(before)
    out["before_auroc"] = auroc(cb, mask[:len(cb)])
    out["before_ns"] = float(before.time[-1] / 1000.0)

    # AFTER: cosolvent
    topo = protein_topology(cos / f"{tag}_cosolv_top.pdb", cos / f"{tag}_protein_top.pdb")
    after = md.load(str(cos / f"{tag}_cosolv.xtc"), top=str(topo))[::STRIDE]
    ca, ser = cavity_p95(after)
    n = min(len(ca), len(mask))
    out["after_auroc"] = auroc(ca[:n], mask[:n])
    out["after_ns"] = float(after.time[-1] / 1000.0)
    out["delta"] = out["after_auroc"] - out["before_auroc"]
    out["mean_cavity_before"] = float(np.mean(cb))
    out["mean_cavity_after"] = float(np.mean(ca))
    return out


if __name__ == "__main__":
    SP = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else
                      "/private/tmp/claude-501/-Volumes-QIDHD-Loki/"
                      "d463e20f-ecd2-46ab-8f02-01cddcf59b6b/scratchpad")
    B = HERE.parent / "benchmark"
    labels = json.loads((B / "cb_labels.json").read_text())
    TARGETS = {"1fd3_A": ("1fd4_G", B / "heldout_cohort/_extracted"),
               "1egw_B": ("6byy_B", B / "heldout_cohort/_extracted"),
               "1u55_A": ("6cww_B", B / "heldout_cohort/_extracted"),
               "2pbk_A": ("1fl1_B", B / "strat_cohort/_extracted")}
    rows = []
    print(f"{'tag':<9} {'lab':>4} {'before':>7} {'after':>7} {'delta':>7} "
          f"{'ns':>6} {'cav_b':>7} {'cav_a':>7}")
    for tag, (le, ex) in TARGETS.items():
        if not (SP / "cosolv_out" / f"{tag}_cosolv.xtc").exists():
            print(f"{tag:<9} (not finished)"); continue
        try:
            r = compare(tag, ex, SP / "cosolv_out", labels, le)
        except Exception as exc:                                  # noqa: BLE001
            print(f"{tag:<9} FAILED {type(exc).__name__}: {exc}"); continue
        rows.append(r)
        print(f"{r['tag']:<9} {r['n_labels']:>4} {r['before_auroc']:>7.3f} "
              f"{r['after_auroc']:>7.3f} {r['delta']:>+7.3f} {r['after_ns']:>6.1f} "
              f"{r['mean_cavity_before']:>7.1f} {r['mean_cavity_after']:>7.1f}")
    if rows:
        (HERE / "cosolvent_result.json").write_text(json.dumps(rows, indent=2))
        d = np.array([r["delta"] for r in rows])
        print(f"\n  mean delta {d.mean():+.3f} over {len(d)} proteins")
