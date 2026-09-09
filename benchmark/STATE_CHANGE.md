# Stage 1: does the method track how much the site actually changes?

The label used everywhere else is binary — a residue either lines a cryptic site
or does not. But crypticity is a *magnitude*: the apo→holo pocket RMSD says how
far the site moved. CryptoBench publishes that as `pRMSD` for all 1,107 entries
(OSF `10.17605/OSF.IO/PZ4A9`), so the better-posed question is whether cavity
score tracks the size of the change rather than a flag.

## What could and could not be tested

**Per-residue state change: not feasible.** Loki's manifest carries
`known_changed_residues` for 222 apo/holo pairs, but only 11 of those apo chains
have an ATLAS trajectory, and after excluding entries already scored and known
holo structures, **2 remain** — one of them with 2 changed residues. There is no
cohort here.

**Per-protein magnitude: feasible**, by joining CryptoBench's `pRMSD` to the
cavity scores already produced. 59 entries, **52 UniProt clusters**, pRMSD
2.01–7.13 Å (median 2.66).

## Result

| | value |
|---|---|
| Spearman rho | **−0.145** (p = 0.31) |
| Pearson r | −0.223 (p = 0.11) |
| mean AUROC, pRMSD below median | 0.649 (n = 26) |
| mean AUROC, pRMSD at/above median | 0.609 (n = 26) |
| difference | −0.040 (p = 0.48) |

**No relationship.** At 52 clusters this has ~80% power for |r| ≥ 0.38, so a
moderate association would have been detected. None was.

## What that means, and it is not flattering

The point estimate is *negative*: the method does slightly worse on the sites
that move most, though not significantly. That direction is mechanistically
expected, and it matters more than its p-value.

A site with a large apo→holo pocket RMSD is precisely one whose apo state looks
least like the bound pocket — the cavity has not formed yet. A site with a small
pRMSD is already most of the way there. So a detector that scores buried void
volume **in the apo ensemble** should find the second kind and miss the first,
which is what the numbers show.

The working conclusion: this measures **partially-formed proto-pockets**, not
conformational crypticity. It ranks sites that are already dented in the apo
structure. The hardest and most valuable cases — sites that genuinely rearrange
to open — are the ones it is weakest on, and no amount of tuning the cavity
detector addresses that, because the signal it needs is absent from the apo
state by definition.

## Status

Exploratory. The join is post-hoc on scores produced for a different purpose,
and `pRMSD` was never used in any tuning or filtering, but this was not
preregistered and is a hypothesis about the method's character rather than a
test of it. The claim that would need prospective testing is the mechanistic
one: that performance falls as apo→holo displacement rises.
