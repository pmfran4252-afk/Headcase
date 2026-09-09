# Preregistration 2 — the ≥10-label stratum, prospectively

Written and pushed before the cohort was downloaded or scored. Nothing below may
change after a result is seen.

## What this tests, and why it needs its own cohort

The first held-out run returned mean AUROC 0.690 overall, and a post-hoc split
showed clusters with ≥10 labelled residues averaging **0.749** against 0.592
below that. That split was chosen after seeing the data, so it is a hypothesis,
not a result. This tests it prospectively.

The mechanism is specific: AUROC estimated from a handful of positives is very
noisy, and CryptoBench entries with few labelled residues are also the ones
whose site annotation is thinnest. If the split is real, the method has a
stated applicability regime. If it is not, the split was noise and should be
dropped rather than quietly retained.

## The cohort

32 entries across **25 UniProt clusters**, built by rules applied before
scoring:

- ≥ 10 labelled residues on the labelled chain (the stratum boundary, fixed here)
- UniProt cluster **never previously scored** — every cluster touched by the
  14-protein tuning set or the 29-entry held-out run is excluded
- ATLAS entry is not a known ligand-bound member of an apo/holo pair
- spatial-clustering join check at z > −2, applied at run time before scoring

The cluster-level exclusion matters and was nearly missed. Candidates such as
`1n8v_B` (labels from `1kx9_A`), `3in9_A` (from `3ikw_A`) and `7e2s_A` (from
`1y6i_A`) are different structures of proteins already in the tuning set. Their
trajectories are unseen but their sites are not, so they are excluded.

## The method

Unchanged from `PREREGISTRATION.md` and untouched since it was frozen: cavity
p95, probe 1.4 Å, enclosure 3, span 8 Å, 3 replicas pooled, stride 5, no
normalisation, no blending. No parameter has been altered in light of the 0.690.

## The prediction

**Mean AUROC ≥ 0.72** across the eligible clusters.

Set below the 0.749 that generated the hypothesis, because a post-hoc slice
regresses — the previous preregistration predicted 0.70 from a tuning figure of
0.827 and observed 0.690, so this one assumes a similar shrink and still asks
for a number above the pooled 0.690.

- **≥ 0.72** — the stratum is a genuine applicability regime; report it as the
  method's stated operating range.
- **0.65 – 0.72** — real but indistinguishable from the pooled 0.690;
  stratification buys nothing and should be dropped.
- **< 0.65** — refuted. The split was noise.

Reported alongside, not optional: SASA and RMSF baselines, the 95% CI, the
cluster count, and the filter outcome.

## What this cannot settle

With 25 clusters and a between-protein SD near 0.15 the standard error is about
0.03, so this can establish that the stratum sits well above chance but **cannot
cleanly separate 0.72 from the pooled 0.690** — those intervals overlap.
Distinguishing them needs roughly 50 clusters, which this data source does not
contain. A result in the 0.65–0.72 band is therefore genuinely inconclusive
about stratification, and will be reported as such rather than argued either way.
