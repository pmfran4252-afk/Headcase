# Cryptic-pocket detection by NMR-style exchange analysis

**Short answer: yes — but not by simulating an NMR experiment.**

The instinct behind the question is right. NMR relaxation dispersion is the best
experimental method there is for finding cryptic pockets, and it works for one
specific reason: it measures *exchange between a visible ground state and a
sparsely populated, structurally invisible excited state*. A cryptic pocket is
exactly such an excited state. So the thing to borrow is the measurement's
logic, not its physics.

Borrowing the physics literally does not work, and it is worth being precise
about why before building anything.

## The timescale gap, which decides the whole design

A CPMG experiment reports exchange broadening `R_ex`, which depends on the
minor population `p_B`, the exchange rate `k_ex`, and the shift difference
`Δω`. You can compute all three from a trajectory. The trouble is the window:

| | window |
|---|---|
| exchange a 100 ns trajectory can **sample** | `k_ex ≳ 10⁸ s⁻¹` |
| exchange a CPMG experiment can **measure** | `k_ex ≲ 10⁵ s⁻¹` |
| overlap | **empty** |

Both numbers are computed, not asserted — `detectability_limit()` gives the
first from the requirement that a fit needs several transitions, and the
Carver–Richards forward model gives the second. Run
`examples/demo_cryptic_sites.py` and you can watch a perfectly well-fitted
excited state (`p_B = 0.29`, `k_ex = 3.6 × 10⁹ s⁻¹`) produce a CPMG dispersion
amplitude of zero, because at that rate the exchange is motionally narrowed to
invisibility.

This has a sharp consequence for scoring. **Do not rank sites on `R_ex`, or on
anything else containing `1/k_ex`.** `k_ex` is the one fitted quantity whose
value cannot be transferred from a short trajectory to the real system, so a
score built on it ranks on a number the simulation got wrong by orders of
magnitude. Rank instead on the **exchange contrast**

```
Φ_ex = p_A · p_B · Δ_obs²
```

which asks how distinct the alternate state is and how much of it there is —
timescale-free, and dominated by the squared amplitude term, which is what
separates a real opening from a loop flickering between two near-identical
states.

The CPMG forward model still earns its place, just in a different role: given a
predicted excited state, `nmr_window()` says which experiment could actually
observe it (CPMG, CEST, R1ρ, HDX-MS, or none). That turns a prediction into a
testable proposal instead of a number.

## Three channels, chosen to fail differently

No single readout survives the timescale gap, so the pipeline runs three and
keeps them visible rather than collapsing them into one number — the *pattern*
of agreement is diagnostic.

| channel | what it measures | blind spot |
|---|---|---|
| **dispersion** | two-state fit → `p_B`, `k_ex`, `Δ` | silent below `k_ex ≈ 10⁸ s⁻¹`, i.e. for most real cryptic sites |
| **tail** | extreme-value extrapolation of pocket volume | extrapolates rather than measures |
| **breathing** | protection-factor anomaly | flags flexible loops alongside real pockets |

**Dispersion** (`exchange.py`) fits a two-state Gaussian HMM per residue and
inverts `P = exp(Q·dt)` exactly for `k_ex`. An HMM rather than a per-frame
mixture, because a mixture assigns frames independently: a noisy frame of state
A that drifts past the midpoint becomes a one-frame excursion into B, and a few
percent of such frames inflate the transition count several-fold. On
ground-truth data the Viterbi path recovers the true transition count *exactly*
(35/35, 9/9, 4/4, 0/0 across the tested range).

**Tail** (`tails.py`) is the answer to the sampling problem. Instead of waiting
for a rare opening, fit a generalised Pareto to the tail of the pocket-volume
distribution and ask what volume corresponds to a state of population `10⁻³`.
Thousands of small excursions in 100 ns constrain how far the large ones reach.
Frames are serially correlated, so the effective sample size comes from
runs-declustering while the tail itself is fit on the marginal — populations
are fractions of *time*, not of independent events. `threshold_stability()`
refits across thresholds; a return level that swings with the threshold is not
evidence.

**Breathing** (`hdx.py`) is the cheapest channel and the only one with a
same-week experimental counterpart. Amide exchange is a rate process, so the
observable is `⟨1/PF⟩`, not `1/⟨PF⟩`, and the gap between the two is precisely
the error you make by reading burial off a single structure. In the test suite
two residues with **identical** mean burial (30.01 vs 29.89 contacts —
indistinguishable in any static or mean-field analysis) separate 128-fold on
this anomaly.

## Does it work?

`examples/demo_cryptic_sites.py` plants three cryptic sites spanning `k_ex =
10⁵` to `2 × 10⁹ s⁻¹`, plus a decoy flexible loop — the false positive this
class of method most often produces.

```
rank  res   score  disp_z  tail_z breath_z chan  truth
   1    8   68.87   32.60  126.24    46.59    3  fast cryptic
   2   17   35.55   -2.38   76.45    38.96    2  slow cryptic
   3   29   28.33    3.55   48.20    40.18    3  glacial cryptic
   …
  39   35                                        flexible loop (decoy)
```

All three real sites rank 1–3, roughly 8× clear of the best false positive, and
the decoy is rejected. Note the second row: the dispersion channel is *negative*
for the slow site, and the tail and breathing channels carry it. That is the
multi-channel design doing the job it exists for.

## Wiring in your models

The pipeline takes arrays, not file formats, so each model plugs in at a
specific point:

- **ATLAS** → the ensemble. Its 3 × 100 ns replicas per protein at 10 ps
  framing are the input to every channel. Pool all three replicas before
  fitting: it triples the trajectory length and thus lowers the dispersion
  channel's detectability floor by 3×.
- **A per-residue embedding model** → `observables.slow_projection()`. Feed a
  `(T, R, D)` embedding trajectory and it returns the slowest direction per
  residue by time-lagged independent component analysis. This is the part that
  is *better* than NMR rather than merely analogous: a chemical shift is only
  accidentally sensitive to cavity formation, whereas a learned coordinate can
  be trained so the contrast `Φ_ex` is maximised for cryptic sites — you are
  designing the spectrometer instead of inheriting one. On synthetic data the
  projection beats the best single raw embedding dimension, and it recovers the
  SNR-limited maximum correlation with the hidden state (0.581 against a
  theoretical 0.59).
- **A per-residue openness or pocket-volume head** → the `openness` argument of
  `score.analyse()`, feeding the tail channel. Any scalar that grows as the site
  opens works: fpocket/mdpocket volume, buried void volume, or a learned score.
- **Structures/coordinates** → `hdx.amide_environment()` for the breathing
  channel.

I could not check how Atlas, Loki and MPSV are represented in this repository —
it was empty at the time of writing — so the adapters are written against those
generic shapes. Point me at the model interfaces and I will wire them directly.

## Validation, before trusting any of this

Synthetic recovery is necessary, not sufficient. The real checks:

1. **CryptoSite** (Cimermancic et al., *JMB* 2016) — ~93 apo/holo pairs where
   the cryptic site is known from the holo structure. Score the apo ensembles
   and measure enrichment of the true site.
2. **PocketMiner's held-out set** (Meller et al., *Nat. Commun.* 2023) — the
   closest published comparator, a GNN predicting cryptic-pocket formation from
   a single structure. Anything proposed here should be measured against it,
   not just against random.
3. **Negative controls from ATLAS itself** — proteins with no known cryptic
   site give the false-positive rate, which is the number that decides whether
   this is usable for triage.
4. **HDX-MS on the top hits** — the breathing channel predicts a measurable
   protection-factor anomaly. This is the cheapest genuine falsification
   available and does not need isotope labelling or a cryoprobe.

## Where this fits

Treat it as **cheap triage, not proof**. The intended architecture is: run all
three channels across ATLAS (cheap, no new simulation), rank, then spend
enhanced sampling — SWISH or mixed-solvent MD with benzene/isopropanol probes,
or metadynamics along a pocket-volume collective variable — only on the top
slice. Mixed-solvent probe occupancy is also the in-silico analogue of
fragment-based NMR screening (STD, WaterLOGSY, ¹⁹F), and it answers the question
none of the three channels here can: whether an openable pocket is a
*druggable* one.

## Honest limitations

- **The dispersion channel will be silent for most real cryptic sites** at
  ATLAS trajectory lengths. This is quantified, not hidden: `detectability_limit()`
  returns the slowest resolvable `k_ex`, so a negative result can be
  distinguished from an absent one.
- **Tail extrapolation to `p = 10⁻⁵` spans three orders of magnitude beyond the
  data.** It orders sites usefully and predicts absolute volumes badly. Use
  `threshold_stability()` and disbelieve unstable fits.
- **Force fields are parameterised on folded ground states** and are expected to
  under-populate open excited states. Every channel inherits that bias.
- **A pocket that opens is not necessarily druggable.** Nothing here scores
  ligandability; that needs the mixed-solvent step.
- **The breathing channel flags flexible loops.** It is the least specific
  channel, which is why it carries the lowest weight and why single-channel hits
  are flagged as such.

## Usage

```python
import numpy as np
from cryptic_dispersion import analyse, project_all, breathing_anomaly

observable, _ = project_all(embeddings)          # (T, R, D) -> (T, R)
ranked = analyse(observable, openness=volumes,   # (T, R)
                 breathing=breathing_anomaly(env), dt_ps=10.0)

for site in ranked[:10]:
    print(site.residue, round(site.score, 2), site.n_channels, site.flags)
```

Requires `numpy` and `scipy` only.

```bash
python tests/test_cryptic_dispersion.py     # 27 tests
python examples/demo_cryptic_sites.py       # end-to-end demonstration
```

## Layout

| file | contents |
|---|---|
| `cryptic_dispersion/exchange.py` | HMM, exact CTMC rate inversion, Carver–Richards & Luz–Meiboom forward models, detectability limit |
| `cryptic_dispersion/tails.py` | declustering, autocorrelation time, GPD tail fit, threshold stability |
| `cryptic_dispersion/hdx.py` | amide environment, Best–Vendruscolo protection factors, breathing anomaly |
| `cryptic_dispersion/observables.py` | per-residue TICA projection, contact numbers, robust scaling |
| `cryptic_dispersion/score.py` | channel z-scoring, consensus ranking, flags |
| `tests/` | 27 tests, each pinned to an analytic limit or a ground-truth generator |
