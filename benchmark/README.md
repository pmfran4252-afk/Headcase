# Real-protein benchmark: CryptoBench labels on ATLAS trajectories

The first time this pipeline meets a protein. Everything before it ran on
generators, and a generator cannot tell you whether a channel measures the thing
it was designed to measure or an artefact of how the data was made.

```sh
./benchmark/fetch_data.sh                    # ~700 MB, not vendored
python benchmark/cryptobench_atlas.py
```

Needs `mdtraj` in addition to numpy and scipy. Results are written after every
entry and an existing result file is resumed, so an interrupted run keeps what
it had.

## Data

**Labels** — CryptoBench (Škrhák et al., *Bioinformatics* 2025), 1,107 apo
structures with cryptic binding-site residues annotated from apo/holo pairs,
keyed by PDB id with residues as `<chain>_<author_resnum>`.

**Trajectories** — ATLAS (Vander Meersche et al., *NAR* 2024), CC-BY-NC. Three
independent 100 ns replicas per chain, frames every 100 ps, 1,001 frames each.

**Intersection: 14 of 1,107.** Found by probing every CryptoBench entry against
the ATLAS API on the chain carrying its labels. That is the practical ceiling on
this benchmark and is worth knowing before designing around it. It is a lower
bound: only the labelled chain was probed.

## Two things that had to be got right

**The three channels read three different measurements.** Burial is a CA contact
count; openness is Shrake–Rupley SASA; breathing uses backbone amide N contacts
and hydrogen bonds. Each is computed independently from the coordinates. They
correlate through the real conformational state, which is the premise of
consensus scoring — but none is an algebraic transform of another, so their
agreement is evidence rather than arithmetic. That is the distinction the
synthetic demo originally failed.

**Replicas are scored separately and averaged, not concatenated.** The project
README recommends pooling all three replicas to triple the trajectory length.
That is right for counting statistics and wrong for every time-series estimate:
concatenating three independent runs inserts two discontinuities that the HMM
reads as transitions and that inflate the integrated autocorrelation time.

## The label join is verified, not assumed

CryptoBench does not document whether its residue numbers are author or CIF
numbering, and a join that lands on the wrong residues produces a confident,
meaningless enrichment. Both are tried, the better match is used, and the match
rate is reported per entry.

The join is then checked by a route that does not involve the numbering at all:
a real binding site is a spatial cluster, so the mean pairwise Cα distance among
mapped residues is compared against random residue sets of the same size drawn
from the same structure.

| entry | mapped | mean pairwise Cα | random | z |
|---|---|---|---|---|
| 1esw_A | 12/12 (author) | 10.3 Å | 27.7 Å | −7.3 |
| 1kx9_A | 15/15 (author) | 11.3 Å | 16.7 Å | −4.4 |
| 2h7g_X | 5/5 (author) | 9.5 Å | 24.4 Å | −3.3 |
| 3b49_A | 20/20 (author) | 16.8 Å | 22.4 Å | −3.7 |

100% of labels map, on author numbering, and every mapped set is significantly
more compact than chance.

## What ATLAS geometry does to the channels

The observability floors are arithmetic, and at 100 ns with 100 ps framing they
are decisive:

| channel | floor at 100 ns, 5% state | verdict |
|---|---|---|
| dispersion | `k_ex ≥ 8.4e8 s⁻¹` | blind — cryptic opening runs near 1e6 |
| tail | `k_ex ≥ 5.3e9 s⁻¹` | blind |
| breathing | population ≥ 0.010 | usable |

So the dispersion channel is gated on its floor *before* the fit rather than
after. Fitting 276–500 Baum–Welch HMMs per replica to populate a channel that
provably cannot speak was the dominant cost of this harness and bought nothing.

This is the honest configuration at ATLAS trajectory lengths, and it is what the
observability work predicted: **breathing is the only channel with reach.**

## Wall status: CryptoBench is spent

All 1,107 CryptoBench apo entries appear in Loki's CryptoBench cohorts (Gate 4,
222 pairs; plus the Gate 6 holdout) — 2,212 PDB entries between them, a 100%
overlap with zero clean entries remaining.

Headcase reads no Loki outcome, so this is not leakage into Loki. The cost is
narrower and real: **results here and Loki's results on these proteins are not
independent evidence and must not later be pooled as though they were.** If an
independent set is needed it has to come from somewhere CryptoBench did not —
CryptoSite's original pairs, PocketMiner's held-out set, or fresh apo/holo pairs
from AHoJ — and the subtraction belongs in the specification, before the freeze.

---

# Result: the pipeline does not beat chance on real proteins

14 proteins, 5,207 residues, 179 labelled cryptic-site residues, one 100 ns
replica each.

```
AUROC : 0.39 0.43 0.44 0.44 0.44 0.48 0.50 0.52 0.52 0.52 0.53 0.60 0.62 0.76
        mean 0.512   median 0.506   sd 0.096   95% CI [0.457, 0.567]
```

| test | result |
|---|---|
| one-sample t-test vs 0.5 | t = +0.47, **p = 0.643** |
| Wilcoxon signed-rank | **p = 1.000** |
| sign test | 7/14 above 0.5, **p = 1.000** |
| median enrichment over base rate | **0.00×** |

**9 of 14 proteins have zero true cryptic residues in the top-k.** The result is
robust to every sensitivity filter — dropping the entry with a bad join, the two
underpowered entries, and the single hit all leave the confidence interval
containing 0.5.

## What was actually tested

Not the three-channel design. It never ran:

| channel | blind on | share |
|---|---|---|
| dispersion | 5,207 / 5,207 residues | **100%** |
| tail | 3,761 / 5,207 residues | **72%** |
| breathing | 0 / 5,207 residues | 0% |

At 100 ns the observability floors rule out dispersion entirely and the tail
channel for most residues, exactly as predicted before the run. So this is a
benchmark of the **protection-factor anomaly alone**, and the finding is that
one channel does not carry the task. Consensus scoring is untested here, because
there was nothing to reach consensus with.

## The one real signal

`2jlq_A` scores AUROC 0.756 with 3.73× enrichment. It survives the correct null:

- Circular-shift permutation of the label mask (preserves site size *and* spatial
  contiguity; residues within a protein are not independent, so Hanley–McNeil is
  the wrong test): null mean 0.501, sd 0.076, **two-sided p = 0.0029**
- Bonferroni over 14 entries: **p = 0.041**
- Its join is the most tightly clustered in the cohort (z = −8.7), so it is not
  a mapping artefact

One protein in fourteen, marginal after correction. Worth a look, not a claim.

## Join verification, all 14

Mean pairwise Cα distance among mapped residues against random residue sets of
the same size. Every entry maps 100% of its labels, but mapping is not
correctness:

| entry | conv | obs Å | rand Å | z |
|---|---|---|---|---|
| 1esw_A | pdb | 10.3 | 27.7 | −7.3 |
| 1hp1_A | pdb | 30.0 | 28.1 | **+0.6 — not clustered** |
| 1kx9_A | pdb | 11.3 | 16.8 | −4.1 |
| 1y6i_A | pdb | 10.6 | 22.6 | −6.6 |
| 2fp1_A | cif | 11.0 | 20.2 | −5.0 |
| 2h7g_X | pdb | 9.5 | 24.5 | −3.5 |
| 2jlq_A | cif | 12.4 | 28.5 | −8.7 |
| 2po4_A | pdb | 11.2 | 39.3 | −10.4 |
| 3b49_A | pdb | 16.8 | 22.5 | −3.5 |
| 3ikw_A | pdb | 14.1 | 25.9 | −4.3 |
| 3vjz_A | pdb | 7.1 | 13.6 | −1.9 (3 labels, underpowered) |
| 4uc8_A | pdb | 12.4 | 15.8 | −0.8 (3 labels, underpowered) |
| 5op0_B | pdb | 16.5 | 24.0 | −4.0 |
| 6irx_A | cif | 11.6 | 30.0 | −5.8 |

`1hp1_A` maps every label and is still not a spatial cluster, which is why the
match rate alone is not sufficient evidence that a join is right.

## What this does and does not establish

It does **not** refute the method. It establishes that at ATLAS trajectory
lengths only one channel can speak, and that channel alone does not find cryptic
sites. Both are consequences of 100 ns being three orders of magnitude short of
where cryptic pockets open.

The informative next step is not more proteins at this length. It is either a
longer-timescale ensemble where the other two channels come off the floor, or
enhanced sampling on a small set — which is what the project README proposed
spending compute on in the first place.
