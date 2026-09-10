# Correction: holo contamination inflated both reported results

## The error

The apo guard excluded ATLAS trajectories that are known ligand-bound members of
an apo/holo pair, using the 222 pairs in Loki's CryptoBench manifest. **The full
CryptoBench dataset names 3,899 holo structures.** The guard was checking a 6%
subset of the list it needed.

A trajectory that is the holo structure starts with the cryptic pocket already
open. Scoring it is not a test of cryptic-site detection, and it scores high for
the wrong reason.

## Effect on the reported numbers

| cohort | as reported | holo-contaminated | **apo-only (corrected)** |
|---|---|---|---|
| held-out | 0.690 (24 clusters) | 0.764 (n=6) | **0.665** (n=18) |
| stratified | 0.640 (14 clusters) | 0.834 (n=3) | **0.588** (n=11) |

Contaminated entries scored 0.10–0.25 higher than clean ones in both cohorts,
consistently and in the expected direction.

**Held-out, corrected:** 0.665, 95% CI [0.559, 0.771], p = 0.0044. The result
stands — still clearly above chance and still ahead of the SASA and RMSF floors.

**Stratified, corrected:** 0.588, 95% CI [0.443, 0.732], p = 0.206. **This no
longer reaches significance.** The stratified cohort does not independently
establish that the method beats chance; it was already REFUTED against its
≥ 0.72 prediction, and the corrected figure weakens it further.

## What does not change

The refutation of the label-count stratification stands and is strengthened —
0.588 is further below the predicted 0.72, not closer to it.

The mechanistic conclusion is reinforced rather than undermined. Trajectories
that begin in the bound state — where the pocket is already formed — are exactly
the ones that scored best. That is the same finding as the pRMSD join and the
skewness result, arriving a third time by a third route: **this method scores
pockets that are already open, and the more open the starting structure, the
better it does.**

## Why it was missed

The guard was written when the only holo list to hand was Loki's manifest, and
the full CryptoBench dataset with all holo ids was not downloaded until the
state-change work, several steps later. The check was never revisited against
the better list once it existed. The lesson is narrow and practical: a guard
built against a partial reference must be re-run when a complete one arrives.

Any future cohort must screen `atlas_entry` against every `holo_pdb_id` in the
CryptoBench dataset, not against a pairing manifest.
