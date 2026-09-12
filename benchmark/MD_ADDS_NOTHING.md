# The molecular dynamics contributes nothing

Same detector, same proteins, same labels. One deposited structure against 100 ns
of MD aggregated as p95:

| | mean AUROC |
|---|---|
| single frame (crystal structure) | **0.707** |
| p95 over the 100 ns ensemble | 0.700 |
| paired difference | **−0.007, p = 0.63** |

Single-frame cavity beats chance on its own (p = 0.0064). The ensemble adds
nothing measurable.

## What this invalidates

Every trajectory-based elaboration in this project. The p95 aggregation, the
extreme-value tail channel, the two-state dispersion fit, the entire
observability apparatus with its detectability floors, replica pooling, the
dynamics observables — none of it beats reading the deposited coordinates once.

It also explains results that looked puzzling in isolation:

- τ, skewness, kurtosis and step asymmetry all failed to add signal
- cavity *level* beat every cavity *change* measure
- the trajectories that scored best were the ones that started in the bound state

There was no dynamical signal to extract. A static geometric descriptor has been
wearing an MD costume.

## What it is worth

The pipeline collapses from about 13 GPU-hours per protein to **0.03 seconds**.
That is a ~1,500,000× reduction and it makes screening the entire PDB feasible.

## The one thing that still needs simulation

Cosolvent MD. Probes cannot bind a static structure, so probe occupancy is the
only observable here that genuinely requires a trajectory — and the only one
whose compute cost is justified.

---

# And the detector fails on real drugged targets

Run on the deposited apo structure of five known drugged or high-value cryptic
sites, using the configuration the control above shows is not a handicap:

| target | AUROC | top-k enrichment |
|---|---|---|
| Bcl-xL BH3 groove | 0.686 | 2.40× |
| ABL1 myristoyl (asciminib, FDA 2021) | 0.611 | 0.00 |
| p53 | 0.569 | 0.00 |
| IL-2 (Ro26-4550) | 0.354 | 0.00 |
| CDK2 allosteric | 0.208 | 0.00 |
| **mean** | **0.486** | **0.48×** |

**At chance**, against 0.707 on the CryptoBench cohort with the identical
detector. Four of five return zero correct residues in the top-k. Asked to find
the pocket that produced asciminib, it does not.

The benchmark flatters the method. Whatever CryptoBench's apo/holo pairs select
for, it is easier than the targets anyone would actually want triaged, and an
AUROC on that cohort should not be read as readiness for target selection.

n = 5, and two of these UniProt ids sit in Loki's spent cohorts, so this is a
demonstration rather than a validated claim. It is still the sharpest available
evidence about whether the tool is ready, and the answer is no.
