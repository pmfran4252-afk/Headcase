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

---

# OUTCOME — recorded after execution, prediction unchanged

**Verdict: PARTIAL.** Predicted mean AUROC ≥ 0.70; observed **0.690**.

29 entries, **24 UniProt clusters** (the unit of analysis — several entries are
the same protein in different crystal forms and are not independent).

| measure | mean | median | 95% CI | >0.5 | p vs 0.5 |
|---|---|---|---|---|---|
| **cavity p95 (frozen)** | **0.690** | 0.745 | [0.601, 0.778] | 19/24 | 1.9e-04 |
| SASA baseline | 0.498 | 0.456 | [0.426, 0.570] | 9/24 | 0.96 |
| RMSF baseline | 0.529 | 0.534 | [0.451, 0.607] | 14/24 | 0.45 |

Paired against its floors: **+0.192 over SASA** (p = 0.0011) and **+0.161 over
RMSF** (p = 0.0060). Wilcoxon against chance p = 4.3e-04; sign test 19/24,
p = 0.0066.

So the effect is real and clearly separated from both baselines, and the
prediction was still missed. Per this document's own terms — "0.60–0.70: real
but weaker than claimed; report the shortfall, do not re-tune and re-report" —
the shortfall is the result.

The drop from the tuning set's 0.827 to 0.690 is the regression this
preregistration anticipated, which is why it predicted 0.70 rather than 0.827.
It anticipated the direction and still overshot the magnitude.

## An applicability signal, logged for a future preregistration only

| stratum | clusters | mean AUROC |
|---|---|---|
| ≥ 10 labelled residues | 15 | 0.749 |
| < 10 labelled residues | 9 | 0.592 |

Spearman against label count rho = +0.35 (p = 0.098), against protein size
rho = +0.40 (p = 0.054).

**This does not rescue the result.** Re-slicing the held-out cohort after seeing
it and reporting 0.749 is precisely the manoeuvre this file exists to prevent.
It is recorded as a hypothesis for a *separate* prospective test with the
stratum boundary fixed in advance, and the headline number stays 0.690.
