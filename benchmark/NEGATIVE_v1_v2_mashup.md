# Negative: combining Headcase cavity with Loki's v2 channels does not help

Recorded so the mashup is not rebuilt in the expectation that it works.

Loki's `cryptic_formula_v2` scores five channels — strain, rewire, burial, face,
lid — from two conformations. An MD trajectory supplies those for free: frame 0
as the unbound state, the most-open frame as the second. Every channel was
computed on the tuning cohort alongside `cavity p95` and combined.

## Individual channels

| channel | verified-7 AUROC |
|---|---|
| **cavity (Headcase)** | **0.802** |
| burial (v2) | 0.716 |
| rewire (v2) | 0.712 |
| strain (v2) | 0.671 |
| lid (v2) | 0.599 |
| face (v2) | 0.563 |

v2's channels are not noise — `rewire` and `burial` carry real signal. Cavity
beats all of them by about 0.09.

## Combinations

| combination | verified-7 |
|---|---|
| cavity + 0.5×rewire | 0.823 |
| **cavity alone** | **0.802** |
| cavity + 0.5×strain | 0.802 |
| cavity + 0.5×face | 0.789 |
| cavity + full v2 | 0.765 |
| cavity + 0.5×lid | 0.732 |
| v2 full sum | 0.698 |

One combination improves on cavity, by **+0.021 on seven proteins**. That is
noise, and it is the same shape as the ≥10-label stratum that read 0.749 and
returned 0.588 under a prospective test.

## The reasoning that failed, which is the useful part

`lid` was expected to complement cavity: it measures a hydrophobic flap
*covering* a void, where cavity measures a void already open. The correlation
data confirms they measure different things — **lid has the lowest correlation
with cavity of any channel, |r| = 0.108**, essentially orthogonal.

Adding it makes performance *worse*, 0.802 → 0.732.

**Low correlation is not complementarity.** Orthogonal information only helps if
it is also informative, and lid at 0.599 is too weak to carry its weight in a
sum. The inference "measures something different, therefore adds value" does not
hold, and it is the inference that motivated the whole mashup.

## What this does not test

v2 expects `predicted_ca` — a *modelled* second conformation, plausibly from a
structure predictor. It was given frame 0 versus the most-open MD frame, because
that is what exists for these proteins. This is a fair substitution but it is not
v2's native input, so these numbers test **"v2's channels computed on MD frames,"
not "v2 as designed."** Its own docstring is correct that it has no validated
performance claim; this does not supply one in either direction.

Cohort note: 9 proteins (7 verified-label) rather than the usual 14/11 — entries
over 400 residues were skipped for the O(n²) distance matrices, and a few for
CA/residue count mismatches. Cavity reads 0.802 here against 0.823 on the full
set, so the subset is not distorting the comparison.
