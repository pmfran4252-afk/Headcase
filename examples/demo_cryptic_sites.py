"""End-to-end demonstration on a synthetic protein with known cryptic sites.

The point is not that the pipeline finds pockets in made-up data.  It is to test
the one claim the architecture rests on: that three channels reading *different*
evidence fail differently, so their agreement carries information.

That claim is only testable if the channels really are separate measurements.
An earlier version of this demo derived pocket volume from the NMR-like
observable and amide burial from the volume, which made all three channels
algebraic transforms of a single latent variable -- for a residue with no
planted site, r(obs, vol) was exactly 1.000.  Agreement was then guaranteed by
construction and measured nothing.

Here each channel is its own instrument: a shared hidden two-state process is
the only thing they have in common, and each reads it through an independent
thermal background and independent measurement noise.  Two further changes make
the test honest:

  * Every residue gets the *same* anharmonic volume response.  Giving the real
    sites an exponential response and the decoy a linear one hands the
    extreme-value channel the answer, since a Pareto fit separates those two
    functional forms whatever they represent.
  * The decoy is a real false positive rather than a weak signal.  A flexible
    surface loop moves *more* than a cryptic site and breathes just as hard; the
    only thing it does not do is open a cavity.  So dispersion and breathing
    both fire on it, and the tail channel alone can reject it.

Three planted sites span the timescale range:

  fast    k_ex ~ 2e9 s^-1   inside the window a 100 ns trajectory can resolve
  slow    k_ex ~ 5e6 s^-1   too slow to fit, but it still partially opens
  glacial k_ex ~ 1e5 s^-1   essentially never opens; only breathing survives
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_dispersion.exchange import detectability_limit, exchange_contrast, fit_two_state
from cryptic_dispersion.hdx import AmideEnvironment, breathing_anomaly
from cryptic_dispersion.score import rank_sites
from cryptic_dispersion.tails import latent_openness, threshold_stability

N_RES, N_FRAMES, DT_PS = 40, 8000, 10.0
TOTAL_NS = N_FRAMES * DT_PS / 1000.0
OBS_NOISE = 0.30
# Thermal breathing of a *closed* pocket, dimensionless.  This has to stay small
# relative to a real opening: the volume response is exponential, so driving it
# with a unit-variance thermal coordinate manufactures a heavy Pareto tail for
# every residue in the protein and blinds the tail channel completely.  A closed
# pocket does not open because a coordinate wandered three sigma.
VOL_THERMAL, VOL_NOISE = 0.35, 0.10

# residue -> (p_B, k_ex, obs amplitude, cavity amplitude, exposure, label)
#
# The decoy moves the most and breathes as hard as any real site.  Its cavity
# amplitude is zero: it flexes without opening anything.
SITES = {
    8:  (0.15, 2.0e9, 2.5, 3.0, 0.85, "fast cryptic"),
    17: (0.05, 5.0e6, 2.5, 3.2, 0.85, "slow cryptic"),
    29: (0.03, 1.0e5, 2.5, 3.5, 0.90, "glacial cryptic"),
    35: (0.25, 3.0e9, 3.0, 0.0, 0.85, "flexible loop (decoy)"),
}


def ou(n: int, tau: float, rng: np.random.Generator) -> np.ndarray:
    """Ornstein-Uhlenbeck process with unit stationary variance."""
    a = np.exp(-1.0 / tau)
    s = np.sqrt(1.0 - a * a)
    x = np.zeros(n)
    e = rng.standard_normal(n)
    for i in range(1, n):
        x[i] = a * x[i - 1] + s * e[i]
    return x


def telegraph(n: int, p_b: float, k_ex: float, rng: np.random.Generator) -> np.ndarray:
    dt_s = DT_PS * 1e-12
    p_ab = 1.0 - np.exp(-p_b * k_ex * dt_s)
    p_ba = 1.0 - np.exp(-(1.0 - p_b) * k_ex * dt_s)
    out = np.zeros(n)
    s = 0
    u = rng.random(n)
    for i in range(n):
        s = (1 if u[i] < p_ab else 0) if s == 0 else (0 if u[i] < p_ba else 1)
        out[i] = s
    return out


def volume_response(u: np.ndarray) -> np.ndarray:
    """Pocket volume from a dimensionless opening coordinate.

    Applied identically to every residue.  A cavity expands super-linearly once
    it starts to open, which is what puts weight in the tail -- but that is a
    property of the readout, not a label, so the decoy and the null residues get
    exactly the same curve.
    """
    return 20.0 + 6.0 * u + 4.0 * np.exp(np.clip(u, None, 4.0)) * (u > 0)


def build(seed: int = 0):
    rng = np.random.default_rng(seed)
    obs = np.empty((N_FRAMES, N_RES))
    vol = np.empty((N_FRAMES, N_RES))
    contacts = np.empty((N_FRAMES, N_RES))
    hbonds = np.empty((N_FRAMES, N_RES))
    hidden = {}

    for i in range(N_RES):
        if i in SITES:
            p_b, k_ex, a_obs, a_vol, a_exp, _ = SITES[i]
            state = telegraph(N_FRAMES, p_b, k_ex, rng)
            hidden[i] = state
        else:
            state = np.zeros(N_FRAMES)
            a_obs = a_vol = a_exp = 0.0

        # Three instruments. Each sees the hidden state through its own thermal
        # background and its own measurement noise, so every correlation between
        # channels has to travel through the state itself.
        obs[:, i] = (a_obs * state
                     + ou(N_FRAMES, rng.uniform(5.0, 30.0), rng)
                     + rng.normal(0.0, OBS_NOISE, N_FRAMES))

        u_vol = (a_vol * state
                 + VOL_THERMAL * ou(N_FRAMES, rng.uniform(5.0, 30.0), rng)
                 + rng.normal(0.0, VOL_NOISE, N_FRAMES))
        vol[:, i] = volume_response(u_vol)

        # Every residue breathes a little, on its own fast background process.
        # Without that heterogeneity the null population is a constant, its MAD
        # collapses, and the median/MAD z-score of any real site runs to four
        # figures -- an artefact of the generator, not a detection.
        base = rng.uniform(0.02, 0.35) * telegraph(
            N_FRAMES, rng.uniform(0.05, 0.25), 1.0e10, rng)
        exposed = np.clip(a_exp * state + base, 0.0, 0.98)
        contacts[:, i] = np.clip(
            rng.normal(30.0, 1.0, N_FRAMES) - 26.0 * exposed
            + 0.8 * ou(N_FRAMES, rng.uniform(5.0, 30.0), rng), 0.0, None)
        hbonds[:, i] = (rng.random(N_FRAMES) > exposed).astype(float)

    return obs, vol, AmideEnvironment(contacts, hbonds), hidden


def independence_check(obs, vol, env) -> None:
    """Show that the channels are separate measurements before trusting agreement."""
    nulls = [i for i in range(N_RES) if i not in SITES]
    def mean_r(a, b, idx):
        return float(np.mean([abs(np.corrcoef(a[:, i], b[:, i])[0, 1]) for i in idx]))
    print("channel independence (mean |r| over the 36 residues with no planted site)")
    print(f"  |r| obs vs volume      {mean_r(obs, vol, nulls):.3f}")
    print(f"  |r| obs vs contacts    {mean_r(obs, env.n_contacts, nulls):.3f}")
    print(f"  |r| volume vs contacts {mean_r(vol, env.n_contacts, nulls):.3f}")
    print("  (all near zero: with no shared hidden state the instruments are unrelated)\n")


def main() -> int:
    obs, vol, env, hidden = build()
    breath = breathing_anomaly(env)
    independence_check(obs, vol, env)

    disp, tails, stable, fits = [], [], [], []
    for i in range(N_RES):
        f = fit_two_state(obs[:, i], dt_ps=DT_PS)
        fits.append(f)
        disp.append(exchange_contrast(f) if f.reliable else float("nan"))
        tails.append(latent_openness(vol[:, i], target_population=1e-3))
        stable.append(bool(threshold_stability(vol[:, i])["stable"]))

    ranked = rank_sites(np.array(disp), np.array(tails), breath,
                        fits=fits, tail_stable=stable)
    order = [s.residue for s in ranked]

    print(f"Synthetic protein: {N_RES} residues, {TOTAL_NS:.0f} ns at {DT_PS:.0f} ps")
    print(f"Dispersion channel can resolve k_ex >= {detectability_limit(TOTAL_NS, 0.05):.1e} s^-1 "
          f"for a 5% state\n")
    print(f"{'rank':>4} {'res':>4} {'score':>7} {'disp_z':>7} {'tail_z':>7} {'breath_z':>8} "
          f"{'chan':>4}  truth")
    print("-" * 78)
    for rank, s in enumerate(ranked[:10], 1):
        truth = SITES.get(s.residue, ("", "", "", "", "", ""))[5] or "-"
        print(f"{rank:>4} {s.residue:>4} {s.score:>7.2f} {s.dispersion_z:>7.2f} "
              f"{s.tail_z:>7.2f} {s.breathing_z:>8.2f} {s.n_channels:>4}  {truth}")

    # A planted site is only findable if it actually opened during the run.
    # Expected transitions = 2 T p_A p_B k_ex; below about one, the site never
    # leaves the closed state and is simply absent from the data.  Reporting a
    # rank for it measures the generator, not the pipeline.
    print("\nplanted site -> observability, then rank")
    print(f"{'res':>5} {'k_ex':>9} {'E[trans]':>9} {'seen':>5} {'rank':>5} {'disp_z':>7} "
          f"{'tail_z':>7} {'breath_z':>8}  label")
    by_res = {s.residue: s for s in ranked}
    t_s = TOTAL_NS * 1e-9
    observable = []
    for res, (p_b, k_ex, _, a_vol, _, label) in sorted(SITES.items()):
        n_exp = 2.0 * t_s * p_b * (1.0 - p_b) * k_ex
        seen = int(np.abs(np.diff(hidden[res])).sum())
        if seen > 0:
            observable.append(res)
        sc = by_res[res]
        print(f"{res:>5} {k_ex:>9.0e} {n_exp:>9.2f} {seen:>5} {order.index(res) + 1:>5} "
              f"{sc.dispersion_z:>7.2f} {sc.tail_z:>7.2f} {sc.breathing_z:>8.2f}  {label}")

    absent = [r for r in SITES if r not in observable]
    if absent:
        print(f"\n  residues {absent} never left the closed state in {TOTAL_NS:.0f} ns.")
        print("  They are absent from the trajectory, so no channel can rank them and")
        print("  they are excluded from the verdict below. Ranking them would score the")
        print("  generator's label, not a detection.")

    real_obs = [r for r in observable if SITES[r][3] > 0]
    decoy_obs = [r for r in observable if SITES[r][3] == 0]
    print("\n  verdict, over the sites that actually occur")
    if not real_obs or not decoy_obs:
        print("    inconclusive: need at least one observable real site and one decoy")
    else:
        worst = max(order.index(r) + 1 for r in real_obs)
        best_decoy = min(order.index(r) + 1 for r in decoy_obs)
        print(f"    worst observable real site  rank {worst}")
        print(f"    best observable decoy       rank {best_decoy}")
        print("    " + ("PASS: every observable real site outranks the decoy"
                        if worst < best_decoy else
                        "FAIL: the decoy outranks an observable real site"))
        margin = by_res[real_obs[0]].score - by_res[decoy_obs[0]].score
        print(f"    margin {margin:+.2f} on a score of {by_res[real_obs[0]].score:.2f}"
              f"  ({100*margin/max(abs(by_res[real_obs[0]].score),1e-9):.0f}%)")

    print("\n  which channel is doing the separating")
    for r in real_obs + decoy_obs:
        sc, lab = by_res[r], SITES[r][5]
        tail_state = "REFUSED (fit failed)" if sc.tail_z == 0.0 else f"{sc.tail_z:+.2f}"
        print(f"    res {r:>2} {lab:<22} dispersion {sc.dispersion_z:>+7.2f}   "
              f"breathing {sc.breathing_z:>+6.2f}   tail {tail_state}")
    print("    The tail channel is the only one that can tell a cavity from a loop,")
    print("    and it carries 0.35 of the weight against 0.65 for the two that cannot.")

    print("\n  sensitivity of the decoy's rank to the channel weights")
    for name, w in (("default     0.40/0.35/0.25", {"dispersion": 0.40, "tail": 0.35, "breathing": 0.25}),
                    ("tail-led    0.20/0.55/0.25", {"dispersion": 0.20, "tail": 0.55, "breathing": 0.25}),
                    ("cavity-only 0.00/1.00/0.00", {"dispersion": 0.00, "tail": 1.00, "breathing": 0.00})):
        o = [x.residue for x in rank_sites(np.array(disp), np.array(tails), breath,
                                           weights=w, fits=fits, tail_stable=stable)]
        if real_obs and decoy_obs:
            wr = max(o.index(r) + 1 for r in real_obs)
            bd = min(o.index(r) + 1 for r in decoy_obs)
            print(f"    {name}   worst real {wr:>2}   decoy {bd:>2}"
                  f"   {'ok' if wr < bd else 'decoy wins'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
