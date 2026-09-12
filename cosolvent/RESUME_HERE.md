# Resume point — occupancy cohort, frozen and not yet run

Everything below is decided and pushed. Nothing has been simulated. Pick up at
"Next action".

## State

`PREREGISTRATION_OCCUPANCY.md` is frozen at commit `17b72b6` and must not be
edited. It fixes:

- **Primary readout:** probe occupancy AUROC (not cavity)
- **Prediction:** mean ≥ 0.60; < 0.55 refutes
- **Secondary:** paired occupancy − cavity > 0, with CI
- **Cohort:** 14 proteins in `benchmark/occupancy_cohort.json`, drawn at random
  from 775 eligible fresh UniProt clusters with `numpy.default_rng(20260911)`
- **Protocol:** CHARMM36m + TIP3P, benzene `BENZ` at 0.25 M, probes placed
  before solvation, PME 1.0 nm, HBonds, hydrogen mass 4 amu, 4 fs, Langevin
  310 K, MC barostat, 100 ps equilibration, 20 ns production

Cohort: `1fvr_A 1x2g_C 1yhv_A 3iec_B 3wb9_C 4j2p_A 4uum_A 4xmq_B 5nia_A 6cqe_B
6dau_E 6ial_F 6qfl_A 7x0f_A`

## Why the cohort could be this large

Earlier cohorts were capped by needing a pre-computed ATLAS trajectory, which
left 28 usable CryptoBench chains in total. Cosolvent MD makes its own
trajectory, so that requirement does not apply and 775 fresh clusters qualify.
GPU time is the constraint now, not data.

## Next action

1. Build the structure-preparation path: fetch from RCSB, strip to the labelled
   chain, drop waters/ligands, add hydrogens. Verified to work on `2xc1_A`
   (9,941 particles under CHARMM36 with no repair needed); `pdbfixer` is
   installed for entries that need residue rebuilding.
2. Verify **one** full system builds and runs a few hundred steps before
   committing GPU to all fourteen. Three topology bugs in this project were
   invisible until analysis, so the analysis path gets exercised first.
3. Launch, daemonized. ~13 h per protein, ~7.5 days total at 20 ns each.
4. Analyse with `cosolvent/compare.py`, which already reports occupancy AUROC,
   cavity, enrichment and the paired difference.

## Do not

- Re-tune anything on the four proteins already run. They are spent.
- Substitute cavity for occupancy as the headline. The preregistration names
  occupancy, and cavity failed its own test (+0.046, p = 0.399).
- Pool any result here with Loki's. All of CryptoBench sits inside Loki's spent
  cohorts.
