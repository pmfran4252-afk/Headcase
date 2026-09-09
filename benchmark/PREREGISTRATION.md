# Preregistration — cavity volume on the held-out cohort

Written before the held-out cohort was scored. Frozen at the commit that adds
this file; nothing below may be changed after seeing a held-out result.

## Status of the 14-protein cohort

**Burned as a tuning set.** Forty-two configurations were evaluated on it:
32 in the probe/aggregation/normalisation sweep, 8 variants of pooling and
consensus blending, and 10 auxiliary-feature combinations. No number measured on
those 14 proteins is evidence for anything any more.

## The frozen method

Per-residue **buried-cavity volume**, `cryptic_dispersion.cavity`:

| parameter | value |
|---|---|
| grid spacing | 1.0 Å |
| probe radius | 1.4 Å (water) |
| `min_enclosure` | 3 (buried on all three axes) |
| `span_A` | 8.0 Å |
| lining cutoff | 5.0 Å |
| frames | 3 replicas pooled, stride 5 |
| aggregation | **p95 across frames** |
| normalisation | none |
| blending | **none — cavity alone** |

Choices that were between-noise on the tuning set and were made on other
grounds, stated so they cannot be relitigated later: p95 (0.827) against max
(0.828) and SD (0.830) is within 0.003; p95 is taken because max is a
single-frame extremum and SD is scale-dependent, while a high quantile is a
return level and matches how the rest of this project talks about tails.

Choices the tuning set actually decided: probe 1.4 Å over 2.0/2.6/3.2 Å
(0.823 vs 0.798/0.726/0.760); no normalisation (uniformly worse); no blending
with breathing (0.827 alone vs 0.819 at 75/25 and 0.790 at 50/50); no auxiliary
feature, since packing density, B-factor, Neq and hydrophobicity each reduced
the verified-set score when added.

## The held-out cohort

Every UniProt-matched candidate in `cohort.json` **excluding the 14 tuning
proteins**, with the pre-specified guards applied:

- exclude any entry whose ATLAS structure is a known ligand-bound member of an
  apo/holo pair (pocket may start open — not a cryptic-detection test)
- exclude any entry whose labels fail the spatial-clustering join check at
  z > −2, using the same test and threshold already applied to the tuning set
- exclude any entry with fewer than 5 labelled residues

The join filter is the one that removed `1hp1_A`, `3vjz_A` and `4uc8_A` from the
tuning set. It was specified before the cavity work and is applied here
unchanged.

## The prediction

**Mean AUROC ≥ 0.70** across the held-out proteins.

Deliberately below the 0.827 measured on the tuning set. Forty-two
configurations were compared there, and the reported figure is the maximum over
that search on eleven proteins, so it is an optimistic estimate of a quantity
that should regress.

- **≥ 0.70** — supported.
- **0.60 – 0.70** — real but weaker than claimed; report the shortfall, do not
  re-tune and re-report.
- **< 0.60** — refuted. The 0.827 was tuning-set overfit, and this file says so
  in advance.

Reported alongside, not optional: mean SASA and RMSF baselines per entry, the
number of proteins, and the smallest effect the cohort size could have detected
at 80% power.

## What this cannot establish

The held-out proteins are still CryptoBench labels on ATLAS trajectories, and
all of CryptoBench sits inside Loki's spent cohorts. A result here is evidence
that cavity volume ranks cryptic-site residues in this data. It is not evidence
of transfer to a different label set, a different force field, or a prospective
target, and it must not be pooled with Loki's results as independent.
