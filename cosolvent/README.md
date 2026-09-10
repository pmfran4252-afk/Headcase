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
