"""Ground-truth generators for validating the exchange analysis."""

from __future__ import annotations

import numpy as np


def two_state_trajectory(
    n_frames: int,
    dt_ps: float,
    p_b: float,
    k_ex: float,
    delta: float = 1.0,
    noise: float = 0.25,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Telegraph process with exactly known ``p_B`` and ``k_ex``.

    Returns ``(observable, state)``.  Rates follow from the two constraints
    ``k_ex = k_AB + k_BA`` and ``p_B = k_AB / k_ex``, so the generator is
    parameterised by precisely the quantities the fitter must recover.
    """
    rng = np.random.default_rng(seed)
    dt_s = dt_ps * 1e-12
    k_ab = p_b * k_ex          # A -> B
    k_ba = (1.0 - p_b) * k_ex  # B -> A
    p_ab = 1.0 - np.exp(-k_ab * dt_s)
    p_ba = 1.0 - np.exp(-k_ba * dt_s)

    state = np.zeros(n_frames, dtype=int)
    s = int(rng.random() < p_b)
    u = rng.random(n_frames)
    for i in range(n_frames):
        if s == 0:
            if u[i] < p_ab:
                s = 1
        elif u[i] < p_ba:
            s = 0
        state[i] = s

    obs = state * delta + rng.normal(0.0, noise, n_frames)
    return obs, state


def embedded_two_state(
    n_frames: int,
    dt_ps: float,
    p_b: float,
    k_ex: float,
    dim: int = 32,
    signal: float = 1.0,
    noise: float = 1.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """The same process hidden in a random rotation of a ``dim``-D embedding.

    This is the realistic case for a learned per-residue representation: the
    exchange is real but no single coordinate carries it.
    """
    rng = np.random.default_rng(seed)
    obs, state = two_state_trajectory(n_frames, dt_ps, p_b, k_ex,
                                      delta=1.0, noise=0.0, seed=seed)
    z = rng.normal(0.0, noise, (n_frames, dim))
    z[:, 0] += signal * obs
    q, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
    return z @ q.T, state


def breathing_site(
    n_frames: int,
    p_open: float,
    closed_contacts: float = 30.0,
    open_contacts: float = 5.0,
    spread: float = 1.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Amide environment for a site whose burial hides rare full openings."""
    rng = np.random.default_rng(seed)
    is_open = rng.random(n_frames) < p_open
    nc = np.where(is_open,
                  rng.normal(open_contacts, spread, n_frames),
                  rng.normal(closed_contacts, spread, n_frames))
    nh = np.where(is_open, 0.0, 1.0)
    return np.clip(nc, 0.0, None), nh


def anharmonic_openness(
    n_frames: int,
    dt_ps: float,
    p_b: float,
    k_ex: float,
    amplitude: float = 3.0,
    tau_ou: float = 15.0,
    seed: int = 0,
) -> np.ndarray:
    """Pocket volume for a soft site that expands super-linearly as it opens.

    Reproduces the situation that defeats a fixed high threshold: the open state
    is well populated *and* long-lived, so the 95th percentile lands far up
    inside its own right tail and only a handful of independent excursions clear
    it, even though the site opens dozens of times.
    """
    rng = np.random.default_rng(seed)
    a = np.exp(-1.0 / tau_ou)
    noise = np.zeros(n_frames)
    e = rng.standard_normal(n_frames)
    for i in range(1, n_frames):
        noise[i] = a * noise[i - 1] + np.sqrt(1.0 - a * a) * e[i]

    _, state = two_state_trajectory(n_frames, dt_ps, p_b, k_ex,
                                    delta=1.0, noise=0.0, seed=seed + 1)
    u = noise + amplitude * state
    return 20.0 + 6.0 * u + 4.0 * np.exp(np.clip(u, None, 4.0)) * (u > 0)
