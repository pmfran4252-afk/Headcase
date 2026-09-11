# Cosolvent MD: make the rare event happen

Everything else in this project computes an observable on unbiased 100 ns MD.
That cannot work for genuinely cryptic sites, and the project measured why three
separate ways:

1. **pRMSD join** — cavity score does not track how far the site actually moves
   (rho = −0.145, n = 52); the point estimate is *negative*.
2. **Dynamics** — cryptic residues show *less* skew, *less* kurtosis and *less*
   exploration than average. They are not the ones making big excursions.
3. **Holo contamination** — the trajectories that scored best were the ones that
   began in the bound state, pocket already open.

All three say the same thing: the detector finds pockets that are already partly
formed, because a cryptic opening at 10³–10⁶ s⁻¹ is simply absent from 100 ns of
unbiased sampling.

A hydrophobic probe does not wait. Benzene partitions into apolar defects and
holds them open, lowering the barrier instead of hoping to sample over it. This
is the mechanism behind CryptoSite's pipeline and behind SWISH.

## No new force-field parameters

CHARMM36 already ships benzene as residue `BENZ` with CGenFF types `CG2R61` /
`HGR61` at ∓0.115 e. An earlier version of `probe.py` defined a duplicate
template from the phenylalanine ring types; OpenMM correctly rejected it as
*"Multiple non-identical matching templates"*. Using the shipped residue means
nothing unvalidated enters the system.

## Protocol

CHARMM36m + TIP3P, PME at 1.0 nm, HBonds constrained, hydrogen mass 4 amu for a
4 fs timestep, Langevin at 310 K, Monte Carlo barostat at 1 bar. Probes placed
before solvation at 0.25 M so water fills around them. 100 ps equilibration,
then production. Only protein atoms are written to the trajectory, so the
existing cavity detector reads the output unchanged.

## The experiment this exists for

Not a cohort — a mechanism test on four proteins that are genuinely cryptic
(high apo→holo pRMSD) and that the method **failed** on:

| trajectory | apo labels | holo | AUROC | pRMSD |
|---|---|---|---|---|
| `1fd3_A` | 1fd4_G | 6cs9 | 0.248 | 7.13 Å |
| `1egw_B` | 6byy_B | 3mu6 | 0.208 | 5.20 Å |
| `1u55_A` | 6cww_B | 4it2 | 0.358 | 4.14 Å |
| `2pbk_A` | 1fl1_B | 4p3h | 0.394 | 3.87 Å |

Run cosolvent MD, run the identical cavity detector, ask one question: **does
the site appear now?** Confirmed if AUROC rises substantially on the same
proteins with the same detector; the diagnosis is wrong if it does not, and the
GPU budget for a full cohort is saved.

Four cases chosen because a mechanism test needs the right cases, not
statistical power.

## Throughput

Measured on an M1 Pro: **68 ns/day** at 36k atoms on the OpenCL platform against
1.1 ns/day on CPU — a 62× difference that decides whether this is feasible at
all. The diagnostic proteins are 41–228 residues, so most run faster.

```sh
python cosolvent/run_cosolvent.py 1fd3_A <extracted_dir> <out_dir> 20 0.25
```

---

# RESULT: the hypothesis is refuted, a different signal appeared

Four proteins, 20 ns benzene at 0.25 M each, ~45 GPU-hours.

```
tag        lab  before   after   delta  occAUC  occ@site  enrich
1fd3_A       7   0.391   0.382  -0.008     n/a       n/a     n/a
1egw_B       7   0.184   0.367  +0.184   0.416     0.499    0.87
1u55_A       8   0.424   0.410  -0.015   0.760     0.553    1.91
2pbk_A      21   0.418   0.440  +0.022   0.702     0.492    2.07
```

## Cosolvent does not rescue cavity detection

Mean delta **+0.046, p = 0.399**. Three of four remain below chance after the
intervention. The proposed mechanism — probes pry the pocket open, the void
detector then sees it — is **not supported**.

## Probe occupancy outperforms cavity on the same trajectories

| | occupancy AUROC | cavity (after) |
|---|---|---|
| 1egw_B | 0.416 | 0.367 |
| 1u55_A | 0.760 | 0.410 |
| 2pbk_A | 0.702 | 0.440 |
| **mean** | **0.626** | **0.406** |

Paired difference **+0.220** (p = 0.133, n = 3 — underpowered, but occupancy
wins in all three). Enrichment at the labelled site reaches 1.91× and 2.07× in
two of three.

`1u55_A` is the clean case: benzene occupied the cryptic pocket nearly twice as
often as the rest of the protein while the pocket never opened enough for a
grid-based detector to register it — cavity 0.410, occupancy 0.760, one
trajectory.

**Probe binding and pocket opening are separable events, and at 20 ns only the
binding is detectable.** The intervention was half right: changing the input
helped, but the quantity being read out was wrong.

## Status

Not a result. n = 3, p = 0.133, and the cases disagree in magnitude (0.760
against 0.416). This is the same shape as the ≥10-label stratum that looked like
0.749 and returned 0.588 when tested properly.

What distinguishes it from a fishing expedition: occupancy was promoted from
control to observable **before any of this data existed**, on the grounds that
FTMap and CryptoSite already rank sites that way. It is a surviving prediction,
not a pattern found by sweeping.

## What a real test looks like

Preregister **probe occupancy** as the primary readout — ≥0.25 M, ≥50 ns,
multiple replicas — on a fresh cohort of UniProt clusters never scored here,
with the threshold stated first. Everything needed exists: the MD setup, the
occupancy metric, the label machinery, the wall and join filters, and the habit
of declaring the number before running.

## A methodological note worth keeping

`2pbk_A` was inspected at 4.8 of 20 ns and reported cavity 0.528 / occupancy
0.612. Finished, it reads 0.440 / 0.702 — both wrong, in opposite directions.
Partial trajectories mislead in whichever direction the early frames happen to
fall, and were flagged unreliable at the time rather than after.
