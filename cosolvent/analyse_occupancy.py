"""Score the preregistered occupancy cohort.

Primary readout is probe occupancy, as fixed in PREREGISTRATION_OCCUPANCY.md.
Cavity p95 is computed on the same trajectory and reported beside it, because
the claim under test is that occupancy beats cavity, and a claim needs its
comparator computed the same way on the same data.

There is no "before" trajectory here. These proteins were selected without an
ATLAS requirement -- cosolvent MD makes its own trajectory, which is what took
the eligible pool from 28 chains to 775 clusters -- so the comparison is
occupancy against cavity within one simulation, not before against after.

Labels are direct: the CryptoBench apo entry and the simulated structure are the
same PDB entry, so author residue numbers map without SIFTS or sequence
alignment. All 14 map at 100% and pass the spatial-clustering join check at
z between -2.1 and -9.1, verified before the run started.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent / "benchmark"))

import mdtraj as md
from cosolvent.compare import cavity_p95, matching_topology, probe_occupancy
from cryptobench_atlas import auroc

STRIDE = 2


def labels_for(pdb_id: str, chain: str, topology, cbd) -> np.ndarray:
    recs = cbd[pdb_id.lower()]
    rec = max(recs, key=lambda r: r["pRMSD"])
    want = {int(s.split("_")[-1]) for s in rec["apo_pocket_selection"]
            if s.split("_")[0] == chain and s.split("_")[-1].lstrip("-").isdigit()}
    idx = {}
    for i, r in enumerate(topology.residues):
        try:
            idx[int(r.resSeq)] = i
        except (ValueError, TypeError):
            pass
    mask = np.zeros(topology.n_residues, bool)
    for w in want:
        if w in idx:
            mask[idx[w]] = True
    return mask, len(want)


def score(tag, out_dir, cbd):
    pid, ch = tag.split("_")
    out = pathlib.Path(out_dir)
    topo = matching_topology(out / f"{tag}_cosolv_top.pdb", out / f"{tag}_cosolv.xtc",
                             out / f"{tag}_protein_top.pdb")
    t = md.load(str(out / f"{tag}_cosolv.xtc"), top=str(topo))[::STRIDE]
    mask, n_want = labels_for(pid, ch, t.topology, cbd)
    if mask.sum() < 3:
        return {"tag": tag, "error": f"only {int(mask.sum())} of {n_want} labels mapped"}
    occ, n_probes = probe_occupancy(t)
    cav, _ = cavity_p95(t)
    m = min(len(occ), len(cav), len(mask))
    at = float(occ[:m][mask[:m]].mean())
    el = float(occ[:m][~mask[:m]].mean())
    return {"tag": tag, "n_residues": int(t.n_residues), "n_labels": int(mask.sum()),
            "n_probes": int(n_probes), "ns": float(t.time[-1] / 1000.0),
            "occupancy_auroc": auroc(occ[:m], mask[:m]),
            "cavity_auroc": auroc(cav[:m], mask[:m]),
            "occ_at_site": at, "occ_elsewhere": el,
            "occ_enrichment": at / max(el, 1e-9)}


def main() -> int:
    SP = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(
        "/private/tmp/claude-501/-Volumes-QIDHD-Loki/"
        "d463e20f-ecd2-46ab-8f02-01cddcf59b6b/scratchpad")
    cbd = json.loads((SP / "cb_dataset.json").read_text())
    cohort = json.loads((HERE.parent / "benchmark" / "occupancy_cohort.json").read_text())
    out_dir = SP / "occ_out"
    rows = []
    print(f"{'entry':<9} {'res':>5} {'lab':>4} {'prb':>4} {'OCC':>7} {'cavity':>7} "
          f"{'enrich':>7} {'ns':>6}")
    print("-" * 60)
    for c in cohort:
        tag = c["apo"]
        if not (out_dir / f"{tag}_cosolv.xtc").exists():
            continue
        try:
            r = score(tag, out_dir, cbd)
        except Exception as exc:                                   # noqa: BLE001
            print(f"{tag:<9} FAILED {type(exc).__name__}: {exc}"); continue
        if "error" in r:
            print(f"{tag:<9} {r['error']}"); continue
        rows.append(r)
        print(f"{r['tag']:<9} {r['n_residues']:>5} {r['n_labels']:>4} {r['n_probes']:>4} "
              f"{r['occupancy_auroc']:>7.3f} {r['cavity_auroc']:>7.3f} "
              f"{r['occ_enrichment']:>7.2f} {r['ns']:>6.1f}")
    if rows:
        (HERE / "occupancy_result.json").write_text(json.dumps(rows, indent=2))
        o = np.array([r["occupancy_auroc"] for r in rows])
        c = np.array([r["cavity_auroc"] for r in rows])
        print(f"\n  n = {len(rows)}   occupancy {o.mean():.3f}   cavity {c.mean():.3f}   "
              f"paired {(o - c).mean():+.3f}")
        print(f"  PREDICTION was >= 0.60; < 0.55 refutes")
        if len(rows) >= 3:
            from scipy import stats
            ci = stats.t.interval(0.95, len(o) - 1, loc=o.mean(), scale=stats.sem(o))
            print(f"  95% CI [{ci[0]:.3f}, {ci[1]:.3f}]  p vs 0.5 = "
                  f"{stats.ttest_1samp(o, 0.5).pvalue:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
