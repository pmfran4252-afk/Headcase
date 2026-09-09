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
