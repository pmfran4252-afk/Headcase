# Cryptic residues are not floppy — they are conformationally ambiguous

Two questions were tested on the tuning cohort. Both are about *dynamics* rather
than shape, and together they say something specific about what a cryptic
residue is.

## 1. Does the character of the jitter carry signal?

`cavity p95` is a high-water mark: the largest the site ever got. It discards
how it got there. A site poised to open should fluctuate slowly and
asymmetrically — big excursions, slow recovery, backtracking — while a rigid one
rattles fast and Gaussian around a fixed shape. Same p95, different dynamics.

Computed on full 1,001-frame trajectories (dynamics need the time resolution
that striding throws away):

| observable | verified-11 AUROC |
|---|---|
| **p95 (incumbent)** | **0.823** |
| p95 × (1+skew) | 0.819 |
| p95 × τ | 0.778 |
| autocorrelation time τ | 0.679 |
| up/down step asymmetry | 0.490 |
| range ÷ sd | 0.440 |
| excess kurtosis | 0.349 |
| skewness | 0.329 |

**Nothing beats magnitude, and combining dynamics with it makes things worse.**

The sub-chance scores are the finding. An AUROC of 0.329 is inverted signal, not
noise: ranking by skewness puts cryptic residues reliably *low*. Same for
kurtosis and range÷sd. Cryptic residues have **less** right-skew, **less** heavy
tails, and explore **less** relative to their own jitter than the average
residue. They are not the ones making dramatic rare excursions.

## 2. Are they well-ordered?

| observable | verified-11 AUROC |
|---|---|
| cavity p95 | 0.823 |
| **pLDDT inverted (low confidence)** | **0.658** |
| B-factor inverted | 0.560 |
| Neq inverted (rigid) | 0.422 |
| pLDDT (well-ordered) | 0.342 |
| cavity + pLDDT | 0.740 |

Cryptic residues have **low** pLDDT, so the literal hypothesis — well-ordered,
confidently predicted — fails. Adding pLDDT to cavity also *hurts*, so the
signal is real but weaker and overlapping.

## The synthesis

Every measure of **thermal motion** puts cryptic residues at or below chance:
RMSF 0.529/0.472, kurtosis 0.349, range÷sd 0.440, skewness 0.329, Neq-rigid
0.422. They are not unusually mobile, on any timescale measured.

The one flexibility-adjacent quantity that carries signal is pLDDT — and **pLDDT
does not measure wobble.** It measures how confidently a single position can be
assigned. Low pLDDT here plausibly reflects *two supported placements* rather
than disorder: the residue is well-ordered in each state and simply is not in
the same place in both. That is the definition of a cryptic site.

So the intuition that these residues are ordered rather than floppy survives;
what fails is equating "ordered" with "high pLDDT". Order on the fast timescale
and ambiguity between states are different things, and only the second is
predictive.

**Caveat that cannot be resolved here:** low pLDDT is also low for genuinely
disordered regions, termini and flexible loops, and cryptic sites often do
involve loop rearrangement. Part of the 0.658 may simply be "cryptic sites live
in loops," which is true but much shallower than the two-state reading. Nothing
in this data separates them.

**The test that would:** AlphaFold with subsampled MSAs. If these residues are
bimodal rather than disordered, reduced-depth AF2 should produce *both*
placements as distinct confident states rather than one smeared uncertain one.

## A bug, recorded

The first pLDDT run read column 0 of `<tag>_pLDDT.tsv`, which is the residue
letter, not the value, and the `float()` failure was swallowed by a bare
`except`. The tell was `pLDDT` and `pLDDT inverted` returning *identical*
AUROCs, which is impossible for a real variable. Silent fallbacks in a numeric
pipeline produce plausible numbers, and a self-consistency check caught it where
the value alone would not have.
