"""Extreme-value extrapolation of pocket opening beyond the simulated window.

The central difficulty in reading cryptic pockets out of ATLAS-style data is a
timescale mismatch: the trajectories are 100 ns, while cryptic openings and the
excited states that NMR relaxation dispersion reports on live at microseconds
to milliseconds.  Waiting for the open state to appear usually fails.

The way around it is to stop treating the open state as an event to be observed
and start treating it as a *tail* to be estimated.  A 100 ns trajectory contains
thousands of small excursions of the pocket volume; extreme-value theory says
the exceedances above a high threshold converge to a generalised Pareto
distribution regardless of the underlying distribution, so those small
excursions constrain how far the large ones reach.

Read the output as a ranking signal, not as a measured population.  A return
level at p = 1e-5 from 1e4 correlated frames is an extrapolation of three orders
of magnitude; it orders sites reliably and predicts absolute volumes poorly.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional, Sequence

import numpy as np
from scipy import stats


def autocorr_time(x: np.ndarray, max_lag: Optional[int] = None) -> float:
    """Integrated autocorrelation time in frames (initial-positive sequence).

    Needed because MD frames are not independent: 10 000 frames at 10 ps carry
    far fewer than 10 000 independent samples, and ignoring that inflates every
    confidence statement made downstream.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 8:
        return 1.0
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var <= 0:
        return 1.0
    max_lag = min(n // 2, max_lag or n // 2)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, nfft)
    acf = np.fft.irfft(f * np.conjugate(f), nfft)[:max_lag].real
    acf /= acf[0]
    tau = 1.0
    for k in range(1, max_lag):
        if acf[k] <= 0.0:
            break
        tau += 2.0 * acf[k]
    return float(max(tau, 1.0))


def decluster(
    x: np.ndarray,
    threshold: float,
    min_gap: int,
) -> tuple[np.ndarray, int]:
    """Runs-declustering: return one maximum per independent excursion.

    Consecutive exceedances of ``threshold`` separated by fewer than ``min_gap``
    sub-threshold frames belong to the same physical opening event and must not
    be counted as independent evidence.
    """
    x = np.asarray(x, dtype=float).ravel()
    above = x > threshold
    if not above.any():
        return np.empty(0), 0

    idx = np.flatnonzero(above)
    splits = np.flatnonzero(np.diff(idx) > min_gap) + 1
    clusters = np.split(idx, splits)
    peaks = np.array([x[c].max() for c in clusters])
    return peaks, len(clusters)


@dataclass
class TailFit:
    threshold: float
    shape: float          # GPD xi; < 0 means a hard upper bound exists
    scale: float          # GPD sigma
    n_clusters: int       # independent excursions -> effective sample size
    n_frames: int
    zeta_u: float         # P(X > threshold), as a fraction of frames
    tau_frames: float     # integrated autocorrelation time
    upper_bound: Optional[float]   # finite only when xi < 0
    ok: bool
    note: str = ""

    def as_dict(self) -> dict:
        return asdict(self)

    def return_level(self, p: float) -> float:
        """Volume exceeded with probability ``p`` (i.e. by a state of population ``p``)."""
        if not self.ok or p <= 0.0 or p >= 1.0:
            return float("nan")
        if p >= self.zeta_u:
            return float(self.threshold)
        ratio = p / self.zeta_u
        if abs(self.shape) < 1e-8:
            return float(self.threshold + self.scale * np.log(1.0 / ratio))
        val = self.threshold + (self.scale / self.shape) * (ratio ** (-self.shape) - 1.0)
        if self.upper_bound is not None:
            val = min(val, self.upper_bound)
        return float(val)

    def exceedance_probability(self, level: float) -> float:
        """Population of a state whose pocket volume reaches ``level``."""
        if not self.ok:
            return float("nan")
        if level <= self.threshold:
            return float(min(1.0, self.zeta_u))
        z = (level - self.threshold) / self.scale
        if abs(self.shape) < 1e-8:
            return float(self.zeta_u * np.exp(-z))
        base = 1.0 + self.shape * z
        if base <= 0.0:
            return 0.0
        return float(self.zeta_u * base ** (-1.0 / self.shape))


def fit_tail(
    x: np.ndarray,
    quantile: float = 0.95,
    min_clusters: int = 25,
    adaptive: bool = True,
) -> TailFit:
    """Fit a generalised Pareto distribution to the upper tail of ``x``.

    ``x`` is a per-frame scalar that grows as the site opens -- pocket volume
    from a grid method, buried-void volume, or a learned openness score.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    tau = autocorr_time(x)
    min_gap = int(max(1, round(tau)))
    u = float(np.quantile(x, quantile))

    _, n_clusters = decluster(x, u, min_gap)
    if adaptive and n_clusters < min_clusters:
        # A well-populated open state pushes the 95th percentile *inside* it,
        # leaving too few excursions above the threshold to fit.  Walk the
        # threshold down until enough independent excursions clear it; this is
        # the common case for a site that opens more than a few percent of the
        # time, and refusing to fit it would silently drop the easiest targets.
        for q in (0.90, 0.85, 0.80, 0.70, 0.60, 0.50):
            if q >= quantile:
                continue
            u_try = float(np.quantile(x, q))
            _, k = decluster(x, u_try, min_gap)
            if k >= min_clusters:
                u, n_clusters = u_try, k
                break
    if n_clusters < min_clusters:
        return TailFit(u, 0.0, 0.0, n_clusters, n, 0.0, tau, None, False,
                       f"only {n_clusters} independent excursions above u")

    # Fit the *marginal* tail on every exceedance, not on cluster maxima.  We
    # want an occupancy -- the fraction of time the site spends open, directly
    # comparable to an NMR p_B -- and cluster maxima are biased high relative to
    # the marginal distribution.  Serial correlation leaves the MLE consistent;
    # it costs effective sample size, which is what `n_clusters` records and
    # what the `min_clusters` gate above enforces.
    excess = x[x > u] - u
    excess = excess[excess > 0]
    if excess.size < min_clusters:
        return TailFit(u, 0.0, 0.0, n_clusters, n, 0.0, tau, None, False,
                       "degenerate excesses")

    shape, _, scale = stats.genpareto.fit(excess, floc=0.0)
    zeta_u = float(np.mean(x > u))
    bound = u - scale / shape if shape < 0 else None

    return TailFit(
        threshold=u,
        shape=float(shape),
        scale=float(scale),
        n_clusters=int(n_clusters),
        n_frames=int(n),
        zeta_u=float(zeta_u),
        tau_frames=float(tau),
        upper_bound=(float(bound) if bound is not None else None),
        ok=True,
    )


def threshold_stability(
    x: np.ndarray,
    quantiles: Sequence[float] = (0.90, 0.925, 0.95, 0.975),
) -> dict:
    """Refit across thresholds; a trustworthy fit is flat in ``xi``.

    Threshold choice is the main free parameter in a peaks-over-threshold
    analysis, and a return level that swings with it is not evidence.  This is
    the diagnostic that decides whether a hit deserves enhanced sampling.
    """
    shapes, levels, used = [], [], []
    for q in quantiles:
        f = fit_tail(x, quantile=q)
        if f.ok:
            shapes.append(f.shape)
            levels.append(f.return_level(1e-3))
            used.append(q)
    if len(shapes) < 2:
        return {"stable": False, "quantiles": used, "shapes": shapes,
                "levels": levels, "spread": float("nan")}
    lv = np.asarray(levels, dtype=float)
    spread = float(np.ptp(lv) / max(abs(np.median(lv)), 1e-12))
    return {
        "stable": bool(spread < 0.35 and np.ptp(np.asarray(shapes)) < 0.4),
        "quantiles": used,
        "shapes": [float(s) for s in shapes],
        "levels": [float(v) for v in lv],
        "spread": spread,
    }


def latent_openness(
    x: np.ndarray,
    target_population: float = 1e-3,
    baseline: str = "median",
) -> float:
    """How much further a site opens at ``target_population`` than it does typically.

    This is the quantity to rank on: it asks "if this site were caught in a
    state as rare as the ones NMR sees, how much bigger would the pocket be?"
    rather than "how big did the pocket get in 100 ns?".
    """
    x = np.asarray(x, dtype=float).ravel()
    fit = fit_tail(x)
    if not fit.ok:
        return float("nan")   # missing, not zero -- see score._z
    base = float(np.median(x)) if baseline == "median" else float(np.mean(x))
    level = fit.return_level(target_population)
    if not np.isfinite(level):
        return float("nan")
    return float(max(level - base, 0.0))
