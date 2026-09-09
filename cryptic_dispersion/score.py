"""Consensus per-residue cryptic-pocket score.

Three channels, deliberately chosen to fail differently:

``dispersion``
    Exchange broadening from a fitted two-state model.  Highest information
    content -- it is the only channel that returns a population and a rate --
    but it requires the excited state to be *visited* within the trajectory, so
    on 100 ns data it is silent for anything slower than about 10^7 s^-1.

``tail``
    Extreme-value extrapolation of the openness observable.  Sees sites that
    only ever partially open, so it survives the undersampling that silences
    the dispersion channel, at the cost of extrapolating rather than measuring.

``breathing``
    Protection-factor anomaly.  Cheapest, most robust, and the only channel
    with a same-week experimental counterpart (HDX-MS).  It reports opening
    without saying how large the resulting pocket is, so it flags flexible
    loops alongside genuine pockets.

A site that scores on all three is worth enhanced sampling.  A site that scores
on one is worth a second look at why the other two are quiet -- the pattern of
disagreement is diagnostic, so :func:`rank_sites` keeps the channels visible
instead of collapsing them into a single number and discarding the reason.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Optional, Sequence

import numpy as np

from .exchange import ExchangeFit, exchange_contrast, fit_two_state
from .observability import ABSENT, SCORED, ResidueStatus, summarise
from .observables import robust_z
from .tails import fit_tail, latent_openness, threshold_stability

DEFAULT_WEIGHTS = {"dispersion": 0.4, "tail": 0.35, "breathing": 0.25}


@dataclass
class SiteScore:
    residue: int
    score: float
    dispersion_z: float
    tail_z: float
    breathing_z: float
    n_channels: int                 # channels agreeing above `channel_cut`
    p_minor: float = float("nan")
    k_ex: float = float("nan")
    reliable_fit: bool = False
    tail_stable: bool = False
    outcome: str = SCORED           # SCORED | NEGATIVE | ABSENT
    flags: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def _z(values: Optional[np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    """Median/MAD z-score plus a per-residue availability mask.

    A channel that could not be computed for a residue is *missing*, not a
    measurement of zero.  Substituting 0.0 before scaling makes a failed fit
    look like a strongly negative result -- and because fits fail most often on
    the sites that open widest, that turns the best candidates into the
    lowest-ranked ones.  Missing entries score 0 *after* scaling, i.e. neutral
    within the z-distribution.

    The mask is returned because neutral-within-the-distribution is not the same
    as neutral-in-the-sum.  A weighted score with a fixed denominator still
    charges the residue that channel's whole weight, which is the identical bug
    one level up; :func:`rank_sites` uses the mask to renormalise per residue.
    """
    if values is None:
        return np.zeros(n), np.zeros(n, dtype=bool)
    v = np.asarray(values, dtype=float)
    good = np.isfinite(v)
    out = np.zeros(n)
    if good.sum() < 3:
        return out, np.zeros(n, dtype=bool)
    out[good] = np.nan_to_num(robust_z(v[good]), nan=0.0, posinf=0.0, neginf=0.0)
    return out, good


def rank_sites(
    dispersion: Optional[np.ndarray] = None,
    tail: Optional[np.ndarray] = None,
    breathing: Optional[np.ndarray] = None,
    residues: Optional[Sequence[int]] = None,
    weights: Optional[dict] = None,
    channel_cut: float = 1.5,
    fits: Optional[Sequence[ExchangeFit]] = None,
    tail_stable: Optional[Sequence[bool]] = None,
    statuses: Optional[Sequence[ResidueStatus]] = None,
) -> list[SiteScore]:
    """Combine per-residue channel values into a ranked list.

    Each channel is converted to a median/MAD z-score before weighting, because
    the channels carry incompatible units (s^-1, cubic angstroms, log units) and
    because most residues in any protein are null -- a median/MAD scale is not
    dragged around by the handful of real hits the way mean/SD is.
    """
    lengths = [len(v) for v in (dispersion, tail, breathing) if v is not None]
    if not lengths:
        raise ValueError("supply at least one channel")
    if len(set(lengths)) != 1:
        raise ValueError(f"channels have inconsistent lengths: {lengths}")
    n = lengths[0]

    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update(weights)
    present = {k: (v is not None) for k, v in
               (("dispersion", dispersion), ("tail", tail), ("breathing", breathing))}

    zd, md = _z(dispersion, n)
    zt, mt = _z(tail, n)
    zb, mb = _z(breathing, n)
    md = md & present["dispersion"]
    mt = mt & present["tail"]
    mb = mb & present["breathing"]

    # Normalise by the weight actually available *for this residue*, not by the
    # weight available for the run.  With a run-level denominator a residue
    # whose dispersion fit failed keeps at most 60% of the score an otherwise
    # identical residue receives -- and dispersion fits fail preferentially on
    # slow exchange, which is what a real cryptic site is.  The penalty was
    # therefore correlated with the label, in the direction that buries the
    # targets.
    avail = w["dispersion"] * md + w["tail"] * mt + w["breathing"] * mb
    total = (w["dispersion"] * zd * md
             + w["tail"] * zt * mt
             + w["breathing"] * zb * mb) / np.where(avail > 0, avail, 1.0)

    residues = list(residues) if residues is not None else list(range(n))
    out: list[SiteScore] = []
    for i in range(n):
        agree = sum(
            int(z[i] > channel_cut)
            for z, ok in ((zd, present["dispersion"]), (zt, present["tail"]),
                          (zb, present["breathing"])) if ok
        )
        flags: list[str] = []
        p_minor = k_ex = float("nan")
        reliable = False
        if fits is not None and i < len(fits) and fits[i] is not None:
            # None is a legitimate entry: it means no fit was attempted, which is
            # what a caller passes when the observability floor rules the
            # dispersion channel blind before paying for a Baum-Welch fit.  That
            # is not the same as a fit that ran and failed, and dereferencing it
            # crashed the first real-protein run on every single entry.
            f = fits[i]
            p_minor, k_ex, reliable = f.p_minor, f.k_ex, f.reliable
            if not f.reliable and f.note:
                flags.append(f"fit: {f.note}")
        outcome = SCORED
        if statuses is not None and i < len(statuses):
            outcome = statuses[i].outcome
            blind = statuses[i].blind_channels
            if blind:
                flags.append("blind: " + ", ".join(blind))
        stable = bool(tail_stable[i]) if tail_stable is not None else False
        if tail_stable is not None and not stable and zt[i] > channel_cut:
            flags.append("tail extrapolation not threshold-stable")
        absent = [name for name, m, ok in
                  (("dispersion", md, present["dispersion"]),
                   ("tail", mt, present["tail"]),
                   ("breathing", mb, present["breathing"]))
                  if ok and not m[i]]
        if absent:
            flags.append("scored without: " + ", ".join(absent))
        if agree == 1 and total[i] > 0:
            flags.append("single-channel hit")

        out.append(SiteScore(
            residue=int(residues[i]), score=float(total[i]),
            dispersion_z=float(zd[i]), tail_z=float(zt[i]), breathing_z=float(zb[i]),
            n_channels=agree, p_minor=float(p_minor), k_ex=float(k_ex),
            reliable_fit=bool(reliable), tail_stable=stable, outcome=outcome,
            flags=flags,
        ))

    # A residue no channel could have seen carries no evidence either way, so it
    # sorts after everything that was actually observed rather than competing on
    # a score that means nothing.  It stays in the list: dropping it would state
    # the ranking conditional on observability, and observability is correlated
    # with the thing being ranked -- the slower a site opens, the more likely it
    # is both a genuine cryptic pocket and invisible at this trajectory length.
    out.sort(key=lambda s: (s.outcome == ABSENT, -s.score))
    return out


def analyse(
    observable: np.ndarray,
    openness: Optional[np.ndarray] = None,
    breathing: Optional[np.ndarray] = None,
    dt_ps: float = 10.0,
    obs_to_rad_s: float = 1.0,
    target_population: float = 1e-3,
    residues: Optional[Sequence[int]] = None,
    weights: Optional[dict] = None,
) -> list[SiteScore]:
    """End-to-end: ``(T, R)`` observable in, ranked sites out.

    ``openness`` is a ``(T, R)`` quantity that grows as the pocket opens (pocket
    volume, buried void volume, a learned openness head).  It defaults to the
    observable itself, which is correct whenever the observable was oriented so
    that larger means more open.
    """
    obs = np.asarray(observable, dtype=float)
    if obs.ndim != 2:
        raise ValueError("observable must be (T, R)")
    _, n_res = obs.shape
    open_arr = obs if openness is None else np.asarray(openness, dtype=float)
    if open_arr.shape != obs.shape:
        raise ValueError("openness must match the observable's shape")

    fits, disp, tail_v, stable = [], [], [], []
    for i in range(n_res):
        f = fit_two_state(obs[:, i], dt_ps=dt_ps, obs_to_rad_s=obs_to_rad_s)
        fits.append(f)
        disp.append(exchange_contrast(f) if f.reliable else float("nan"))
        tail_v.append(latent_openness(open_arr[:, i], target_population))
        stable.append(bool(threshold_stability(open_arr[:, i])["stable"]))

    return rank_sites(
        dispersion=np.array(disp), tail=np.array(tail_v),
        breathing=None if breathing is None else np.asarray(breathing, float),
        residues=residues, weights=weights, fits=fits, tail_stable=stable,
    )
