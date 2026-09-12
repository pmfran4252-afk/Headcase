# Preregistration 3 — probe occupancy as the primary readout

Written and pushed before any structure was downloaded or any simulation run.
Nothing below may change after a result is seen.

## What is being tested

On four proteins, cosolvent MD failed to improve cavity detection (mean delta
+0.046, p = 0.399) but **probe occupancy outperformed cavity on the same
trajectories**: 0.626 against 0.406, winning in all three cases with occupancy
data, with enrichment at the labelled site of 1.91× and 2.07×.

`1u55_A` is the mechanism in one protein: benzene occupied the cryptic pocket
nearly twice as often as the rest of the surface while the pocket never opened
enough for a void detector to register it — occupancy 0.760, cavity 0.410, one
trajectory. **Probe binding and pocket opening are separable events, and at
these timescales only the binding is detectable.**

That is n = 3 and the cases disagree (0.760 against 0.416). This tests it
properly.

## The cohort, and why it is no longer tiny

Every previous cohort was bounded by a requirement that has now been dropped.
CryptoBench entries needed a *pre-computed ATLAS trajectory*, which left 28
usable chains in the entire dataset and, after exclusions, nothing fresh.
Cosolvent MD generates its own trajectory, so the ATLAS requirement does not
apply and **775 fresh UniProt clusters** become eligible.

**14 proteins drawn at random** from those 775 with `numpy.default_rng(20260911)`
— seed fixed here, selection made before any structure was inspected. No
filtering on size, difficulty, pRMSD or label count beyond the eligibility rules
below, so the draw cannot be steered.

Eligibility, all applied before the draw:

- ≥ 10 labelled residues on the chain
- UniProt cluster never scored in any previous cohort here
- the apo entry is not itself a holo structure anywhere in CryptoBench (checked
  against all 3,899 holo ids, not the 222-pair subset that caused an earlier
  contamination error)

Drawn cohort: `1fvr_A 1x2g_C 1yhv_A 3iec_B 3wb9_C 4j2p_A 4uum_A 4xmq_B 5nia_A
6cqe_B 6dau_E 6ial_F 6qfl_A 7x0f_A` — pRMSD median 2.85 (range 2.02–6.10),
labels median 25.

## Method, frozen

Structures are fetched from RCSB, stripped to the labelled chain, hydrogens
added, then CHARMM36m + TIP3P with **benzene (residue `BENZ`) at 0.25 M**,
probes placed before solvation, PME at 1.0 nm, HBonds constrained, hydrogen mass
4 amu, 4 fs step, Langevin 310 K, Monte Carlo barostat, 100 ps equilibration,
**20 ns production**. Identical to the four already run.

**Primary readout: probe occupancy** — per-residue fraction of frames with a
benzene carbon within 5.0 Å, AUROC against the labels. Reported alongside:
cavity p95 on the same trajectory, occupancy enrichment at the labelled site,
and the paired occupancy − cavity difference.

## The prediction

**Mean occupancy AUROC ≥ 0.60** across the cohort.

Set below the 0.626 seen at n = 3 because small-sample estimates in this project
have regressed consistently: 0.827 → 0.690 and 0.749 → 0.640, drops of 0.11–0.14.
A prediction of 0.60 already assumes most of that shrink.

- **≥ 0.65** — strongly supported; occupancy is the readout this method should use.
- **0.55 – 0.65** — supported but modest; better than cavity, not a solved problem.
- **< 0.55** — refuted. The n = 3 result was noise and cosolvent occupancy is
  not the answer.

Secondary, and the claim that actually matters: **paired occupancy − cavity > 0**
on the same trajectories, reported with its confidence interval.

## Stated in advance: what this cannot settle

At n = 14 with a between-protein SD near 0.15 the standard error is ≈ 0.040, so
the 95% interval is about ±0.086. This can separate 0.65 from chance. It **cannot**
separate 0.60 from 0.65, and it cannot establish an effect smaller than ≈ 0.13
above chance at 80% power.

20 ns is short. It is defensible for occupancy specifically — probes equilibrate
in nanoseconds, unlike pocket opening — but a null at this length is evidence
about *this protocol*, not about cosolvent methods generally, and will be
reported that way.

Labels remain CryptoBench, which sits entirely inside Loki's spent cohorts, so
nothing here may be pooled with Loki results as independent evidence.
