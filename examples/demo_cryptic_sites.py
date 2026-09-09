"""End-to-end demonstration on a synthetic protein with known cryptic sites.

The point of the demo is not that the pipeline finds pockets in made-up data --
it is to show *which channel* finds *which kind* of site, because that is the
design question.  Three planted sites span the timescale range:

  fast    k_ex ~ 2e9 s^-1   inside the window a 100 ns trajectory can resolve
  slow    k_ex ~ 5e6 s^-1   too slow to fit, but it still partially opens
  glacial k_ex ~ 1e5 s^-1   essentially never opens; only breathing survives

plus a decoy: a flexible surface loop that breathes without forming a pocket,
which is the false positive this class of method most often produces.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_dispersion.exchange import detectability_limit, exchange_contrast, fit_two_state, nmr_window
from cryptic_dispersion.hdx import AmideEnvironment, breathing_anomaly
from cryptic_dispersion.score import rank_sites
from cryptic_dispersion.tails import latent_openness, threshold_stability

N_RES, N_FRAMES, DT_PS = 40, 8000, 10.0
TOTAL_NS = N_FRAMES * DT_PS / 1000.0

# residue -> (p_B, k_ex, opening amplitude, anharmonicity, label)
SITES = {
    8:  (0.15, 2.0e9, 3.0, 0.9, "fast cryptic"),
    17: (0.05, 5.0e6, 3.5, 1.0, "slow cryptic"),
    29: (0.03, 1.0e5, 4.0, 1.1, "glacial cryptic"),
    35: (0.25, 3.0e9, 0.4, 0.0, "flexible loop (decoy)"),
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


def build(seed: int = 0):
    rng = np.random.default_rng(seed)
    obs = np.empty((N_FRAMES, N_RES))
    vol = np.empty((N_FRAMES, N_RES))
    contacts = np.empty((N_FRAMES, N_RES))
    hbonds = np.empty((N_FRAMES, N_RES))

    for i in range(N_RES):
        u = ou(N_FRAMES, tau=rng.uniform(5.0, 30.0), rng=rng)
        amp, anh = 0.0, 0.0
        if i in SITES:
            p_b, k_ex, amp, anh, _ = SITES[i]
            u = u + amp * telegraph(N_FRAMES, p_b, k_ex, rng)

        obs[:, i] = u
        # An anharmonic (soft) site expands super-linearly once it starts to
        # open, which is what puts weight in the tail even when the fully open
        # state is never reached.  A rigid site responds linearly.
        vol[:, i] = 20.0 + 6.0 * u + anh * 4.0 * np.exp(np.clip(u, None, 4.0)) * (u > 0)
        opening = np.clip((vol[:, i] - 20.0) / 12.0, 0.0, None)
        contacts[:, i] = np.clip(rng.normal(30.0, 1.0, N_FRAMES) - 7.0 * opening, 0.0, None)
        hbonds[:, i] = (opening < 0.8).astype(float)

    return obs, vol, AmideEnvironment(contacts, hbonds)


def main() -> int:
    obs, vol, env = build()
    breath = breathing_anomaly(env)

    disp, tails, stable, fits = [], [], [], []
    for i in range(N_RES):
        f = fit_two_state(obs[:, i], dt_ps=DT_PS)
        fits.append(f)
        disp.append(exchange_contrast(f) if f.reliable else float("nan"))
        tails.append(latent_openness(vol[:, i], target_population=1e-3))
        stable.append(bool(threshold_stability(vol[:, i])["stable"]))

    ranked = rank_sites(np.array(disp), np.array(tails), breath,
                        fits=fits, tail_stable=stable)

    print(f"Synthetic protein: {N_RES} residues, {TOTAL_NS:.0f} ns at {DT_PS:.0f} ps")
    print(f"Dispersion channel can resolve k_ex >= {detectability_limit(TOTAL_NS, 0.05):.1e} s^-1 "
          f"for a 5% state\n")
    print(f"{'rank':>4} {'res':>4} {'score':>7} {'disp_z':>7} {'tail_z':>7} {'breath_z':>8} "
          f"{'chan':>4}  truth")
    print("-" * 78)
    for rank, s in enumerate(ranked[:10], 1):
        truth = SITES.get(s.residue, ("", "", "", "", ""))[4] or "-"
        print(f"{rank:>4} {s.residue:>4} {s.score:>7.2f} {s.dispersion_z:>7.2f} "
              f"{s.tail_z:>7.2f} {s.breathing_z:>8.2f} {s.n_channels:>4}  {truth}")

    order = [s.residue for s in ranked]
    print("\nplanted site -> recovered rank")
    for res, (p_b, k_ex, _, _, label) in sorted(SITES.items()):
        f = fits[res]
        lim = detectability_limit(TOTAL_NS, p_b)
        seen = "yes" if k_ex >= lim else "no (too slow)"
        print(f"  res {res:>3} {label:<22} k_ex={k_ex:.0e}  rank {order.index(res) + 1:>2}"
              f"   dispersion resolvable: {seen}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
