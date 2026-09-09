"""Turning model output into an NMR-like observable.

A CPMG experiment is sensitive to cryptic pockets only by accident: a chemical
shift responds to ring currents and packing, and a cavity opening under an
aromatic happens to move it.  That accident is why real dispersion data is hard
to interpret -- a large ``|dw|`` tells you something moved, not what.

In silico we are not stuck with the accident.  Any per-residue, per-frame
scalar can play the role of the shift, so we can choose one that is sensitive
to the thing we care about.  Given a learned per-residue embedding (the natural
output of a structure model such as Loki or MPSV), the best scalar is the
*slowest* direction in that embedding: conformational exchange is by
construction the slow process, so the coordinate with the longest
autocorrelation is the one carrying the exchange signal, and everything
orthogonal to it is fast thermal noise that would only broaden the states.

:func:`slow_projection` is a one-component time-lagged independent component
analysis (TICA) run per residue, which solves exactly that problem.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Projection:
    """One residue's synthetic shift trajectory and its provenance."""

    series: np.ndarray        # (T,) unit-variance observable
    weights: np.ndarray       # (D,) direction in embedding space
    eigenvalue: float         # lag-tau autocorrelation of the projection
    timescale_frames: float   # implied timescale, -lag / ln(eigenvalue)
    rank: int                 # embedding dimensions kept after whitening


def slow_projection(
    z: np.ndarray,
    lag: int = 10,
    var_cutoff: float = 0.99,
    ridge: float = 1e-8,
) -> Projection:
    """Project a per-residue embedding trajectory onto its slowest direction.

    Parameters
    ----------
    z
        ``(T, D)`` embedding trajectory for one residue.
    lag
        Lag in frames.  It should exceed the vibrational decorrelation time but
        stay far below the exchange time; for 10 ps ATLAS frames, 10-50 frames
        (100-500 ps) is a sensible window.
    """
    z = np.asarray(z, dtype=float)
    if z.ndim != 2:
        raise ValueError("z must be (T, D)")
    t, d = z.shape
    if t <= lag + 2:
        raise ValueError("trajectory shorter than the lag")

    zc = z - z.mean(axis=0, keepdims=True)
    c0 = (zc.T @ zc) / t
    a, b = zc[:-lag], zc[lag:]
    ctau = (a.T @ b + b.T @ a) / (2.0 * (t - lag))

    # Whiten on C(0); embeddings are typically rank-deficient, so drop the
    # near-null directions rather than inverting through them.
    c0 = c0 + ridge * np.trace(c0) / max(d, 1) * np.eye(d)
    evals, evecs = np.linalg.eigh(c0)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    keep = evals > 0
    if var_cutoff < 1.0 and keep.any():
        frac = np.cumsum(evals[keep]) / np.sum(evals[keep])
        n_keep = int(np.searchsorted(frac, var_cutoff) + 1)
        keep = np.zeros(d, dtype=bool)
        keep[:n_keep] = True
    rank = int(keep.sum())
    if rank == 0:
        raise ValueError("embedding has no variance")

    w_mat = evecs[:, keep] / np.sqrt(evals[keep])
    m = w_mat.T @ ctau @ w_mat
    m = 0.5 * (m + m.T)
    lam, vec = np.linalg.eigh(m)
    top = int(np.argmax(lam))
    direction = w_mat @ vec[:, top]

    series = zc @ direction
    sd = float(series.std())
    if sd > 0:
        series = series / sd
        direction = direction / sd

    ev = float(np.clip(lam[top], -0.999999, 0.999999))
    ts = float(-lag / np.log(abs(ev))) if 0.0 < abs(ev) < 1.0 else float("inf")
    return Projection(series, direction, ev, ts, rank)


def project_all(
    z: np.ndarray,
    lag: int = 10,
    var_cutoff: float = 0.99,
) -> tuple[np.ndarray, list[Projection]]:
    """Apply :func:`slow_projection` to every residue of a ``(T, R, D)`` tensor.

    Returns the ``(T, R)`` observable matrix and the per-residue projections.
    """
    z = np.asarray(z, dtype=float)
    if z.ndim != 3:
        raise ValueError("z must be (T, R, D)")
    t, r, _ = z.shape
    out = np.empty((t, r), dtype=float)
    projections: list[Projection] = []
    for i in range(r):
        p = slow_projection(z[:, i, :], lag=lag, var_cutoff=var_cutoff)
        out[:, i] = p.series
        projections.append(p)
    return out, projections


def contact_number(
    coords: np.ndarray,
    resid: np.ndarray,
    cutoff: float = 6.5,
    exclude_window: int = 2,
    chunk: int = 256,
) -> np.ndarray:
    """Per-residue heavy-atom contact number, a physical fallback observable.

    Use this when no learned embedding is available: burial responds directly to
    a cavity opening, so it works, just with less contrast than a trained
    coordinate.

    ``coords`` is ``(T, A, 3)`` in nanometres or angstroms -- ``cutoff`` must
    match.  ``resid`` is the ``(A,)`` residue index of each atom.
    """
    coords = np.asarray(coords, dtype=float)
    resid = np.asarray(resid, dtype=int)
    t, a, _ = coords.shape
    residues = np.unique(resid)
    out = np.zeros((t, residues.size), dtype=float)
    c2 = cutoff * cutoff

    for k, r in enumerate(residues):
        sel = resid == r
        near = np.abs(resid - r) <= exclude_window
        partners = ~near
        if not partners.any():
            continue
        for start in range(0, t, chunk):
            stop = min(start + chunk, t)
            centre = coords[start:stop][:, sel, :].mean(axis=1)
            d2 = ((coords[start:stop][:, partners, :] - centre[:, None, :]) ** 2).sum(-1)
            out[start:stop, k] = (d2 < c2).sum(axis=1)
    return out


def standardize(x: np.ndarray, axis: int = 0) -> np.ndarray:
    """Zero-mean, unit-variance along ``axis``; constant columns pass through."""
    x = np.asarray(x, dtype=float)
    mu = x.mean(axis=axis, keepdims=True)
    sd = x.std(axis=axis, keepdims=True)
    return (x - mu) / np.where(sd > 0, sd, 1.0)


def robust_z(x: np.ndarray) -> np.ndarray:
    """Median/MAD z-score -- the right scale for ranking when most sites are null."""
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    scale = 1.4826 * mad
    if scale <= 0:
        scale = x.std() or 1.0
    return (x - med) / scale
