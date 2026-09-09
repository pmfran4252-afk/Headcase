"""Validation suite.

Every physics assertion here was checked against an analytic limit or a
ground-truth generator before being written down; these tests exist to keep it
that way.  Runnable with pytest or directly (``python tests/test_cryptic_dispersion.py``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_dispersion.exchange import (  # noqa: E402
    detectability_limit, dispersion_amplitude, exchange_contrast, fit_hmm2,
    fit_two_state, nmr_window, r2eff_carver_richards, r2eff_luz_meiboom,
    rates_from_transition_matrix, viterbi2,
)
from cryptic_dispersion.hdx import (  # noqa: E402
    AmideEnvironment, breathing_anomaly, ensemble_ln_pf, rate_averaged_ln_pf,
)
from cryptic_dispersion.observables import robust_z, slow_projection  # noqa: E402
from cryptic_dispersion.score import rank_sites  # noqa: E402
from cryptic_dispersion.tails import autocorr_time, decluster, fit_tail  # noqa: E402
from tests.synthetic import (  # noqa: E402
    anharmonic_openness, breathing_site, embedded_two_state, two_state_trajectory,
)

NU = np.array([25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1000.0])


# ---------------------------------------------------------------- dispersion

def test_carver_richards_matches_luz_meiboom_in_fast_exchange():
    """The exact two-site solution must reduce to the fast-exchange limit."""
    for p_b in (0.5, 0.10, 0.03):
        cr = r2eff_carver_richards(NU, p_b, 5.0e4, 1.0e3)
        lm = r2eff_luz_meiboom(NU, p_b, 5.0e4, 1.0e3)
        scale = max(float(lm.max()), 1e-30)
        assert np.max(np.abs(cr - lm)) / scale < 1e-3, f"p_B={p_b}"


def test_dispersion_is_quenched_at_high_pulsing_rate():
    """Infinitely fast refocusing removes exchange broadening entirely."""
    for p_b, k_ex, dw in [(0.03, 1e3, 3e3), (0.5, 5e4, 1e3), (0.01, 2e2, 5e3)]:
        assert abs(r2eff_carver_richards(np.array([1e7]), p_b, k_ex, dw)[0]) < 1e-3


def test_slow_exchange_plateau_is_p_b_times_k_ex():
    """In slow exchange R_ex saturates at the escape rate, not at phi_ex/k_ex."""
    p_b, k_ex, dw = 0.02, 200.0, 2.0e4
    plateau = r2eff_carver_richards(np.array([1e-6]), p_b, k_ex, dw)[0]
    assert abs(plateau - p_b * k_ex) / (p_b * k_ex) < 1e-3


def test_dispersion_profile_is_monotone_and_finite():
    """cosh(eta+) overflows float64 at low nu; the log-domain branch must hold."""
    prof = r2eff_carver_richards(np.logspace(-2, 5, 400), 0.03, 1500.0, 3000.0)
    assert np.all(np.isfinite(prof))
    assert np.all(np.diff(prof) <= 1e-9)


def test_dispersion_amplitude_vanishes_without_exchange():
    flat = fit_two_state(np.zeros(2000) + 1e-9, dt_ps=10.0)
    assert abs(dispersion_amplitude(flat)) < 1e-6


# ------------------------------------------------------------------ exchange

def test_viterbi_recovers_the_true_transition_count():
    """State assignment must be exact, or every rate downstream is wrong."""
    for p_b, k_ex in [(0.30, 1e9), (0.10, 1e9), (0.03, 1e9), (0.30, 1e8)]:
        obs, state = two_state_trajectory(10000, 10.0, p_b, k_ex,
                                          delta=1.0, noise=0.25, seed=7)
        h = fit_hmm2(obs)
        path = viterbi2(obs, h["mu"], h["sigma"], h["trans"])
        truth = int(np.count_nonzero(np.diff(state)))
        got = int(np.count_nonzero(np.diff(path)))
        assert got == truth, f"p_B={p_b} k_ex={k_ex:.0e}: {got} vs {truth}"


def test_recovers_ground_truth_when_exchange_is_resolvable():
    """Within the detectable window, p_B and k_ex come back correct."""
    obs, _ = two_state_trajectory(20000, 10.0, 0.20, 1e9, delta=1.0,
                                  noise=0.25, seed=3)
    f = fit_two_state(obs, dt_ps=10.0)
    assert f.reliable, f.note
    assert abs(f.p_minor - 0.20) / 0.20 < 0.25
    assert abs(f.k_ex - 1e9) / 1e9 < 0.35


def test_undersampled_exchange_is_flagged_not_reported():
    """A silent channel must announce itself; this is the whole safety story."""
    obs, _ = two_state_trajectory(10000, 10.0, 0.03, 1e7, delta=1.0,
                                  noise=0.25, seed=7)
    f = fit_two_state(obs, dt_ps=10.0)
    assert not f.reliable and f.note


def test_pure_noise_does_not_produce_a_confident_two_state_fit():
    rng = np.random.default_rng(11)
    f = fit_two_state(rng.standard_normal(8000), dt_ps=10.0)
    assert not f.reliable


def test_ctmc_inversion_round_trips():
    """P = exp(Q dt) inversion must return the rate that generated it."""
    for k_ex in (1e7, 1e8, 1e9, 1e10):
        for p_b in (0.05, 0.2, 0.5):
            dt = 1e-11
            a = p_b * (1.0 - np.exp(-k_ex * dt))
            b = (1.0 - p_b) * (1.0 - np.exp(-k_ex * dt))
            got, ok = rates_from_transition_matrix(np.array([[1 - a, a], [b, 1 - b]]), dt)
            assert ok and abs(got - k_ex) / k_ex < 1e-6


def test_detectability_limit_scales_correctly():
    """Longer trajectories resolve slower exchange, linearly."""
    assert detectability_limit(300.0, 0.05) < detectability_limit(100.0, 0.05)
    assert np.isclose(detectability_limit(100.0, 0.05) / detectability_limit(1000.0, 0.05),
                      10.0, rtol=1e-9)
    assert detectability_limit(100.0, 0.01) > detectability_limit(100.0, 0.30)


# ---------------------------------------------------------------------- tail

def test_gpd_tail_reproduces_the_empirical_exceedance_curve():
    """The extrapolation is only worth anything if it is calibrated where we can check."""
    rng = np.random.default_rng(1)
    n, phi = 40000, 0.9
    e = rng.standard_normal(n)
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    x = np.abs(x)
    f = fit_tail(x)
    assert f.ok
    for level in (5.0, 6.0, 7.0):
        emp = float((x > level).mean())
        if emp > 5e-4:
            assert abs(f.exceedance_probability(level) - emp) / emp < 0.45, level


def test_return_level_increases_as_population_decreases():
    rng = np.random.default_rng(5)
    x = np.abs(rng.standard_normal(20000))
    f = fit_tail(x)
    assert f.return_level(1e-2) < f.return_level(1e-4) < f.return_level(1e-6)


def test_autocorrelation_time_detects_serial_dependence():
    rng = np.random.default_rng(2)
    assert autocorr_time(rng.standard_normal(20000)) < 2.0
    n, phi = 20000, 0.9
    e = rng.standard_normal(n)
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    assert autocorr_time(x) > 5.0


def test_declustering_collapses_one_excursion_to_one_event():
    x = np.zeros(1000)
    x[100:130] = 5.0      # one long excursion
    x[500:505] = 5.0      # a second, well separated
    _, n = decluster(x, threshold=1.0, min_gap=10)
    assert n == 2


def test_tail_fit_refuses_when_there_are_too_few_excursions():
    assert not fit_tail(np.concatenate([np.zeros(999), [10.0]])).ok


# ----------------------------------------------------------------------- hdx

def test_breathing_anomaly_separates_sites_with_identical_mean_burial():
    """The discriminating case: static burial is equal, dynamics are not."""
    nc_static, nh_static = breathing_site(6000, 0.0, closed_contacts=30.0, seed=1)
    nc_breath, nh_breath = breathing_site(6000, 0.03, closed_contacts=30.8, seed=2)
    env = AmideEnvironment(np.stack([nc_static, nc_breath], 1),
                           np.stack([nh_static, nh_breath], 1))
    mean_burial = env.n_contacts.mean(axis=0)
    assert abs(mean_burial[0] - mean_burial[1]) < 0.5      # indistinguishable statically
    anomaly = breathing_anomaly(env)
    assert anomaly[1] > 10.0 * max(anomaly[0], 1e-6)       # separated dynamically


def test_rate_averaging_never_exceeds_mean_field_protection():
    """1/<exp(-lnPF)> <= <lnPF> by Jensen; a violation means a sign error."""
    rng = np.random.default_rng(4)
    env = AmideEnvironment(rng.uniform(0, 35, (500, 12)),
                           (rng.random((500, 12)) < 0.6).astype(float))
    assert np.all(rate_averaged_ln_pf(env) <= ensemble_ln_pf(env) + 1e-9)


# --------------------------------------------------------------- observables

def test_slow_projection_beats_the_best_raw_embedding_dimension():
    z, state = embedded_two_state(8000, 10.0, 0.15, 2e9, dim=32,
                                  signal=1.0, noise=1.0, seed=6)
    proj = slow_projection(z, lag=20)
    got = abs(np.corrcoef(proj.series, state)[0, 1])
    best = max(abs(np.corrcoef(z[:, j], state)[0, 1]) for j in range(z.shape[1]))
    assert got > best


def test_robust_z_is_not_dragged_by_a_few_large_hits():
    x = np.concatenate([np.zeros(200), [50.0, 60.0]])
    assert abs(np.median(robust_z(x))) < 1e-9


# --------------------------------------------------------------------- score

def test_rank_sites_puts_the_multi_channel_hit_first():
    n = 40
    disp = np.zeros(n); tail = np.zeros(n); breath = np.zeros(n)
    disp[7] = tail[7] = breath[7] = 10.0     # agrees on all three
    tail[19] = 12.0                          # louder, but alone
    ranked = rank_sites(disp, tail, breath)
    assert ranked[0].residue == 7
    assert ranked[0].n_channels == 3
    assert any("single-channel" in f for f in
               next(s for s in ranked if s.residue == 19).flags)


def test_rank_sites_rejects_mismatched_channel_lengths():
    try:
        rank_sites(np.zeros(10), np.zeros(11))
    except ValueError:
        return
    raise AssertionError("inconsistent channel lengths must raise")


def test_exchange_contrast_ranks_amplitude_over_flexibility():
    """A loop flickering between near-identical states must not outrank an opening."""
    opening, _ = two_state_trajectory(6000, 10.0, 0.08, 2e9, delta=3.0,
                                      noise=0.25, seed=2)
    loop, _ = two_state_trajectory(6000, 10.0, 0.45, 2e9, delta=0.35,
                                   noise=0.25, seed=2)
    assert exchange_contrast(fit_two_state(opening, dt_ps=10.0)) > \
        5.0 * exchange_contrast(fit_two_state(loop, dt_ps=10.0))


def test_contrast_is_independent_of_the_fitted_rate():
    """Unlike R_ex, the ranking statistic must not carry MD's wrong k_ex."""
    fast, _ = two_state_trajectory(8000, 10.0, 0.20, 2e9, delta=2.0, noise=0.2, seed=1)
    slow, _ = two_state_trajectory(8000, 10.0, 0.20, 4e8, delta=2.0, noise=0.2, seed=1)
    ff, fs = fit_two_state(fast, dt_ps=10.0), fit_two_state(slow, dt_ps=10.0)
    assert ff.reliable and fs.reliable
    assert ff.k_ex > 2.0 * fs.k_ex                       # rates genuinely differ
    a, b = exchange_contrast(ff), exchange_contrast(fs)
    assert abs(a - b) / max(a, b) < 0.25                 # contrast does not


def test_nmr_window_reports_md_timescales_as_unreachable():
    """MD-rate exchange is motionally narrowed; saying otherwise misleads a wet lab."""
    assert "none" in nmr_window(3.6e9)
    assert "CPMG" in nmr_window(1e4)
    assert dispersion_amplitude(fit_two_state(
        two_state_trajectory(6000, 10.0, 0.2, 2e9, delta=3.0, noise=0.2, seed=1)[0],
        dt_ps=10.0)) < 1e-3


def test_failed_channel_scores_neutral_not_negative():
    """A channel that could not be computed must not sink the residue."""
    n = 20
    tail = np.full(n, 5.0)
    tail[3] = np.nan                       # failed fit
    ranked = rank_sites(tail=tail)
    missing = next(s for s in ranked if s.residue == 3)
    assert missing.tail_z == 0.0
    assert all(np.isfinite(s.score) for s in ranked)


def test_adaptive_threshold_fits_a_widely_populated_open_state():
    """A site that opens often and widely must not be dropped for opening too well.

    With a 15% open state reached in long episodes, the 95th percentile sits far
    up inside the open state's own tail and only a handful of independent
    excursions clear it.  Refusing to fit there would discard exactly the sites
    the pipeline exists to find, so the threshold walks down until the data
    supports a fit.
    """
    vol = anharmonic_openness(8000, 10.0, p_b=0.15, k_ex=3.9e9, seed=4)
    strict = fit_tail(vol, adaptive=False)
    relaxed = fit_tail(vol, adaptive=True)
    assert not strict.ok, "expected the fixed 95th-percentile fit to fail here"
    assert relaxed.ok and relaxed.n_clusters > strict.n_clusters
    assert relaxed.threshold < strict.threshold


def _main() -> int:
    tests = [(k, v) for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())
