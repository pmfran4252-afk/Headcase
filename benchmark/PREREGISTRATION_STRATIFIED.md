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

---

## Amendment 1 — label transfer by alignment (made before any score)

The cohort as first built could not be scored: 28 of 32 entries mapped **zero**
labels. The cause is a data gap, not the method. SIFTS `pdb_chain_uniprot`
leaves `PDB_BEG` empty for most rows and the PDBe API returns
`author_residue_number: None` for the same entries, so no author-to-UniProt
offset exists for those labelled chains. Tracing it out:

```
CryptoBench chains with a SIFTS UniProt  : 607
  ...and >=10 labels                     : 399
  ...and any ATLAS entry for that UniProt:  28
  of those 27 clusters, already spent    :  26
  fresh remaining via SIFTS              :   1   (and its ATLAS entry is known holo)
```

Via SIFTS the source is exhausted. The 25 fresh clusters are real proteins with
real trajectories whose labels simply cannot be placed by that route.

**Change:** labels are transferred by aligning the label chain to the ATLAS
chain directly (`align_labels.py`), which needs no UniProt intermediary. Two
crystal forms of one protein align at near-identity, so the correspondence is
exact where it exists. Two guards, fixed here:

- **sequence identity ≥ 0.90** between the two chains, else the pairing is
  rejected as not the same protein
- the spatial-clustering join check at z > −2 is unchanged and still applied
  afterwards; it remains independent of how the mapping was produced, and on the
  five entries used to validate the aligner it still rejected one (`3bl9_A`,
  z = −0.1) while the other four came back between −6.6 and −9.7

**Unchanged:** the scoring method, every cohort rule, and the prediction of
mean AUROC ≥ 0.72.

**Why this is not outcome-driven:** no AUROC has been computed for any entry in
this cohort. The run was stopped during the filter pass, before the scoring loop
was reached, and the only numbers seen were mapped-label counts and join
z-scores — both inputs to eligibility, neither a result. A fix that restores
data the method never saw cannot select for a favourable answer.

---

# OUTCOME — recorded after execution, prediction unchanged

**Verdict: REFUTED.** Predicted mean AUROC ≥ 0.72; observed **0.640**, which
falls in the band this document defined in advance as "the split was noise."

16 entries, **14 UniProt clusters** — the preregistration assumed 25, and half
the cohort was lost to the join check and identity guard.

| measure | mean | median | 95% CI | >0.5 | p vs 0.5 |
|---|---|---|---|---|---|
| **cavity p95 (frozen)** | **0.640** | 0.683 | [0.514, 0.767] | 10/14 | 0.033 |
| SASA baseline | 0.479 | 0.467 | [0.401, 0.556] | 6/14 | 0.56 |
| RMSF baseline | 0.480 | 0.461 | [0.377, 0.584] | 6/14 | 0.69 |

Paired: **+0.162 over SASA** (p = 0.037), **+0.160 over RMSF** (p = 0.054).

## What is refuted, and what is not

**Refuted: stratifying by label count.** The ≥10-label stratum scored 0.749 on
the held-out data, and this test predicted ≥ 0.72 allowing for regression. It
came back at 0.640 — *below* the pooled 0.690, not above it. There is no
evidence that label count identifies an applicability regime, and the 0.749
should be treated as a slice of noise. Reporting it as the method's operating
range would have been wrong, which is what a prospective test is for.

**Not refuted: the method.** Cavity remains above chance (p = 0.033) and ahead
of both baselines by ~0.16, on fourteen clusters that were never used for
anything else. That is consistent with the held-out 0.690 and with the method
having a single performance level rather than two regimes.

## The honest caveat, and its limit

The interval [0.514, 0.767] spans the refutation band, the pooled 0.690 and the
0.72 threshold, so this cohort cannot cleanly place the stratum anywhere. That
was stated before the run, when the expected cohort was 25 clusters; at 14 it is
worse.

That is a reason to distrust the *precision* of 0.640, not a reason to avoid the
verdict. The prediction was a point threshold, the point estimate missed it in
the unfavourable direction, and the pre-specified band is the band. Treating a
wide interval as grounds to withhold a refutation would make the preregistration
unfalsifiable, which is the opposite of why it exists.

## One inconsistency worth recording

The stratum boundary was applied to labels **available on the source chain**,
but scoring uses labels **successfully mapped onto the trajectory**. One cluster
entered with ≥10 source labels and only 8 mapped. A future version should apply
the boundary post-mapping.
