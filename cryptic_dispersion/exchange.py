"""Two-state conformational exchange analysis on MD-derived observables.

This is the computational transposition of an NMR relaxation-dispersion
experiment.  A CPMG experiment reports, per residue, on exchange between a
visible ground state A and a sparsely populated "invisible" excited state B.
The information content is three numbers -- the minor population ``p_B``, the
exchange rate ``k_ex`` and the shift difference ``|dw|`` -- from which the
exchange broadening ``R_ex`` follows.

Cryptic pockets are exactly such excited states, so the same three numbers are
what we want out of a trajectory.  The difference (and the point of doing this
in silico) is that the observable does not have to be a chemical shift, and the
structure of state B is not invisible: we know which frames are in it.

Conventions
-----------
``dw``      angular frequency difference, rad/s
``k_ex``    s^-1
``nu_cpmg`` Hz
``tau_cp``  1 / (4 * nu_cpmg), s
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np


# --------------------------------------------------------------------------
# 1-D two-component Gaussian mixture (EM), numpy only.
# --------------------------------------------------------------------------

def _gauss(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    sigma = max(float(sigma), 1e-12)
    z = (x - mu) / sigma
    return np.exp(-0.5 * z * z) / (sigma * np.sqrt(2.0 * np.pi))


def fit_gmm2(
    x: np.ndarray,
    n_iter: int = 300,
    tol: float = 1e-8,
    seed: int = 0,
) -> dict:
    """Fit a two-component 1-D Gaussian mixture by expectation-maximisation.

    Returns a dict with ``mu``, ``sigma``, ``weight`` (each length 2, ordered by
    increasing ``mu``) and the per-frame posterior ``resp`` of shape ``(T, 2)``.
    """
    x = np.asarray(x, dtype=float).ravel()
    if x.size < 4:
        raise ValueError("need at least 4 samples to fit a two-state model")

    rng = np.random.default_rng(seed)
    lo, hi = np.percentile(x, [15.0, 85.0])
    if hi - lo < 1e-12:
        hi = lo + (np.std(x) + 1e-9)
    mu = np.array([lo, hi], dtype=float)
    sigma = np.full(2, max(np.std(x), 1e-9) / 2.0)
    w = np.array([0.5, 0.5])

    prev = -np.inf
    for _ in range(n_iter):
        dens = np.stack([w[k] * _gauss(x, mu[k], sigma[k]) for k in range(2)], axis=1)
        tot = dens.sum(axis=1)
        # A component can collapse; nudge it back rather than dividing by zero.
        bad = tot <= 1e-300
        if bad.any():
            dens[bad] = 0.5
            tot[bad] = 1.0
        resp = dens / tot[:, None]

        nk = resp.sum(axis=0) + 1e-12
        w = nk / x.size
        mu = (resp * x[:, None]).sum(axis=0) / nk
        var = (resp * (x[:, None] - mu) ** 2).sum(axis=0) / nk
        sigma = np.sqrt(np.maximum(var, 1e-18))
        # Guard against a zero-width component locking onto a single frame.
        floor = 1e-3 * (np.std(x) + 1e-12)
        sigma = np.maximum(sigma, floor)

        ll = float(np.log(np.maximum(tot, 1e-300)).sum())
        if abs(ll - prev) < tol * max(1.0, abs(ll)):
            break
        prev = ll

    order = np.argsort(mu)
    mu, sigma, w = mu[order], sigma[order], w[order]
    dens = np.stack([w[k] * _gauss(x, mu[k], sigma[k]) for k in range(2)], axis=1)
    tot = np.maximum(dens.sum(axis=1), 1e-300)
    resp = dens / tot[:, None]
    _ = rng  # seed reserved for future multi-start initialisation
    return {"mu": mu, "sigma": sigma, "weight": w, "resp": resp}


# --------------------------------------------------------------------------
# Core-set (milestoning) state assignment.
# --------------------------------------------------------------------------

def core_assign(resp: np.ndarray, core_cut: float = 0.85) -> np.ndarray:
    """Assign frames to states using cores, to avoid counting recrossings.

    A frame belongs to a core when its posterior exceeds ``core_cut``.  Frames
    in the transition region inherit the last core visited.  Counting
    transitions between cores instead of between raw maximum-posterior labels
    is what keeps ``k_ex`` from being inflated by barrier-top jitter -- the
    single most common way an MD rate estimate goes wrong.
    """
    resp = np.asarray(resp, dtype=float)
    labels = np.full(resp.shape[0], -1, dtype=int)
    labels[resp[:, 0] > core_cut] = 0
    labels[resp[:, 1] > core_cut] = 1

    seen = np.where(labels >= 0)[0]
    if seen.size == 0:
        return np.full(resp.shape[0], -1, dtype=int)

    out = labels.copy()
    last = labels[seen[0]]
    for i in range(seen[0], out.size):
        if out[i] < 0:
            out[i] = last
        else:
            last = out[i]
    out[: seen[0]] = labels[seen[0]]
    return out


# --------------------------------------------------------------------------
# Two-state Gaussian hidden Markov model.
# --------------------------------------------------------------------------

def fit_hmm2(
    x: np.ndarray,
    n_iter: int = 80,
    tol: float = 1e-9,
    seed: int = 0,
) -> dict:
    """Baum-Welch fit of a two-state Gaussian HMM to one observable trajectory.

    An HMM rather than a per-frame mixture, because the quantity we most need to
    get right is the *rate*, and a mixture assigns each frame independently: a
    frame of state A whose noise carries it past the midpoint becomes a
    one-frame excursion into B, and a few percent of such frames inflate the
    transition count several-fold.  The HMM's transition prior makes a brief
    excursion expensive to explain, so it is attributed to noise instead --
    which is what it is.
    """
    x = np.asarray(x, dtype=float).ravel()
    t = x.size
    if t < 20:
        raise ValueError("need at least 20 frames to fit an HMM")

    init = fit_gmm2(x, seed=seed)
    mu = init["mu"].copy()
    sigma = np.maximum(init["sigma"].copy(), 1e-9)
    pi = np.clip(init["weight"].copy(), 1e-6, None)
    pi /= pi.sum()
    trans = np.array([[0.99, 0.01], [0.01, 0.99]])

    floor = 1e-3 * (np.std(x) + 1e-12)
    prev = -np.inf
    gamma = np.zeros((t, 2))
    xi_sum = np.zeros((2, 2))

    for _ in range(n_iter):
        b = np.stack([_gauss(x, mu[k], sigma[k]) for k in range(2)], axis=1)
        b = np.maximum(b, 1e-300)

        # Forward pass with per-frame scaling (avoids underflow over 1e4 frames).
        alpha = np.empty((t, 2))
        scale = np.empty(t)
        a0 = pi * b[0]
        scale[0] = a0.sum()
        alpha[0] = a0 / scale[0]
        for i in range(1, t):
            ai = (alpha[i - 1] @ trans) * b[i]
            scale[i] = ai.sum()
            if scale[i] <= 0:
                scale[i] = 1e-300
            alpha[i] = ai / scale[i]

        beta = np.empty((t, 2))
        beta[-1] = 1.0
        for i in range(t - 2, -1, -1):
            beta[i] = (trans @ (b[i + 1] * beta[i + 1])) / scale[i + 1]

        gamma = alpha * beta
        gamma /= np.maximum(gamma.sum(axis=1, keepdims=True), 1e-300)

        xi_sum[:] = 0.0
        for k in range(2):
            for l in range(2):
                xi_sum[k, l] = np.sum(
                    alpha[:-1, k] * trans[k, l] * b[1:, l] * beta[1:, l] / scale[1:]
                )

        trans = xi_sum / np.maximum(xi_sum.sum(axis=1, keepdims=True), 1e-300)
        nk = gamma.sum(axis=0) + 1e-12
        mu = (gamma * x[:, None]).sum(axis=0) / nk
        var = (gamma * (x[:, None] - mu) ** 2).sum(axis=0) / nk
        sigma = np.maximum(np.sqrt(np.maximum(var, 1e-18)), floor)
        pi = np.clip(gamma[0], 1e-12, None)
        pi /= pi.sum()

        ll = float(np.log(scale).sum())
        if abs(ll - prev) < tol * max(1.0, abs(ll)):
            break
        prev = ll

    order = np.argsort(mu)
    mu, sigma = mu[order], sigma[order]
    trans = trans[np.ix_(order, order)]
    gamma = gamma[:, order]

    # Stationary distribution of the 2x2 chain, in closed form.
    a, b_off = trans[0, 1], trans[1, 0]
    denom = a + b_off
    pops = np.array([b_off / denom, a / denom]) if denom > 0 else np.array([0.5, 0.5])

    return {"mu": mu, "sigma": sigma, "trans": trans, "pops": pops,
            "gamma": gamma, "loglik": prev}


def viterbi2(x: np.ndarray, mu: np.ndarray, sigma: np.ndarray,
             trans: np.ndarray) -> np.ndarray:
    """Most likely state path (log-domain), used only for transition counting."""
    x = np.asarray(x, dtype=float).ravel()
    t = x.size
    log_b = np.stack([np.log(np.maximum(_gauss(x, mu[k], sigma[k]), 1e-300))
                      for k in range(2)], axis=1)
    log_a = np.log(np.maximum(trans, 1e-300))

    delta = np.full((t, 2), -np.inf)
    psi = np.zeros((t, 2), dtype=int)
    delta[0] = np.log(0.5) + log_b[0]
    for i in range(1, t):
        for l in range(2):
            cand = delta[i - 1] + log_a[:, l]
            psi[i, l] = int(np.argmax(cand))
            delta[i, l] = cand[psi[i, l]] + log_b[i, l]

    path = np.zeros(t, dtype=int)
    path[-1] = int(np.argmax(delta[-1]))
    for i in range(t - 2, -1, -1):
        path[i] = psi[i + 1, path[i + 1]]
    return path


def rates_from_transition_matrix(trans: np.ndarray, dt_s: float) -> tuple[float, bool]:
    """Exact continuous-time ``k_ex`` from a discretely sampled 2x2 chain.

    The frame-to-frame matrix is ``P = exp(Q dt)``, whose non-unit eigenvalue is
    ``1 - a - b = exp(-k_ex dt)``.  Inverting that -- rather than dividing a
    transition count by a dwell time -- is what keeps the estimate unbiased when
    the frame spacing is not negligible against the dwell times, which on 10 ps
    ATLAS frames it often is not.
    """
    a, b = float(trans[0, 1]), float(trans[1, 0])
    lam = 1.0 - a - b
    if lam <= 0.0 or lam >= 1.0:
        # lam <= 0 means the chain flips faster than it is sampled (aliased);
        # lam >= 1 means no transitions were seen at all.
        return 0.0, False
    return float(-np.log(lam) / dt_s), True


# --------------------------------------------------------------------------
# Result container.
# --------------------------------------------------------------------------

@dataclass
class ExchangeFit:
    """Per-residue outcome of the in-silico dispersion analysis."""

    p_minor: float          # population of the sparsely populated state
    k_ex: float             # s^-1
    delta_obs: float        # (minor - major) in the units of the observable
    dw: float               # rad/s, after unit conversion
    phi_ex: float           # p_A p_B dw^2, rad^2/s^2
    r_ex: float             # phi_ex / k_ex, s^-1  (fast-exchange limit)
    n_transitions: int
    separation: float       # state overlap diagnostic (Cohen's d)
    minor_is_upper: bool    # True when the minor state has the larger mean
    reliable: bool
    note: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def fit_two_state(
    x: np.ndarray,
    dt_ps: float,
    obs_to_rad_s: float = 1.0,
    min_transitions: int = 8,
    min_separation: float = 1.2,
    seed: int = 0,
) -> ExchangeFit:
    """Extract ``p_B``, ``k_ex`` and ``dw`` from one observable trajectory.

    Parameters
    ----------
    x
        Observable time series for a single residue, one value per frame.
    dt_ps
        Frame spacing in picoseconds (ATLAS ships 10 ps).
    obs_to_rad_s
        Multiplier converting the observable's units to rad/s.  For a predicted
        1H shift at 600 MHz this is ``2*pi*600e6*1e-6`` rad/s per ppm; leave at
        1.0 to work in the observable's own units, which is the right choice for
        a learned latent coordinate where "ppm" has no meaning.
    min_transitions, min_separation
        Reliability gates.  Both failure modes are real and common: a 100 ns
        trajectory cannot resolve a rate from two crossings, and a two-state fit
        to a unimodal distribution will report a confident, entirely fictitious
        minor state.  Neither is silently returned -- ``reliable`` is set False
        and ``note`` says which gate failed.
    """
    x = np.asarray(x, dtype=float).ravel()
    dt_s = float(dt_ps) * 1e-12
    if dt_s <= 0:
        raise ValueError("dt_ps must be positive")

    hmm = fit_hmm2(x, seed=seed)
    mu, sigma, trans, pops = hmm["mu"], hmm["sigma"], hmm["trans"], hmm["pops"]

    pooled = np.sqrt(0.5 * (sigma[0] ** 2 + sigma[1] ** 2))
    separation = float(abs(mu[1] - mu[0]) / max(pooled, 1e-12))

    k_ex, rate_ok = rates_from_transition_matrix(trans, dt_s)

    path = viterbi2(x, mu, sigma, trans)
    n_trans = int(np.count_nonzero(np.diff(path)))

    minor = int(np.argmin(pops))
    major = 1 - minor
    p_minor = float(pops[minor])
    delta_obs = float(mu[minor] - mu[major])
    dw = float(abs(delta_obs) * obs_to_rad_s)
    phi_ex = float(pops[0] * pops[1] * dw * dw)
    r_ex = float(phi_ex / k_ex) if k_ex > 0 else 0.0

    notes = []
    if not rate_ok:
        notes.append("no resolvable exchange in the sampled window")
    if n_trans < min_transitions:
        notes.append(f"only {n_trans} state transitions")
    if separation < min_separation:
        notes.append(f"states overlap (d={separation:.2f})")
    if p_minor * x.size < 5:
        notes.append("minor state occupies fewer than 5 frames")

    return ExchangeFit(
        p_minor=p_minor,
        k_ex=k_ex,
        delta_obs=delta_obs,
        dw=dw,
        phi_ex=phi_ex,
        r_ex=r_ex,
        n_transitions=n_trans,
        separation=separation,
        minor_is_upper=bool(minor == 1),
        reliable=not notes,
        note="; ".join(notes),
    )


# --------------------------------------------------------------------------
# Forward models: what a spectrometer would have measured.
# --------------------------------------------------------------------------

def r2eff_luz_meiboom(
    nu_cpmg: np.ndarray,
    p_b: float,
    k_ex: float,
    dw: float,
    r2_0: float = 0.0,
) -> np.ndarray:
    """Fast-exchange CPMG dispersion profile (Luz-Meiboom).

    ``R2eff = R2_0 + (phi_ex/k_ex) * [1 - (4 nu / k_ex) tanh(k_ex / 4 nu)]``
    """
    nu = np.asarray(nu_cpmg, dtype=float)
    if k_ex <= 0:
        return np.full_like(nu, r2_0)
    phi = p_b * (1.0 - p_b) * dw * dw
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = k_ex / (4.0 * nu)
        damp = np.where(nu > 0, (4.0 * nu / k_ex) * np.tanh(arg), 0.0)
    return r2_0 + (phi / k_ex) * (1.0 - damp)


def _arccosh_cr(
    d_plus: float,
    d_minus: float,
    eta_plus: np.ndarray,
    eta_minus: np.ndarray,
) -> np.ndarray:
    """``arccosh(D+ cosh(eta+) - D- cos(eta-))``, evaluated without overflowing.

    ``eta_plus`` grows as ``1/nu``, so at the low-field end of a dispersion
    profile ``cosh(eta_plus)`` overflows float64 while the function itself stays
    perfectly finite.  Above the crossover we therefore evaluate the logarithm
    directly, using ``arccosh(S) -> ln(2S)`` for large ``S``:

        ln(2S) = eta+ + ln(D+) + ln1p(e^-2eta+ - 2 (D-/D+) cos(eta-) e^-eta+)

    ``D+ >= 1/2`` always holds (its numerator ``psi + 2 dw^2`` equals
    ``xi^2 + dw^2 + 4 p_A p_B k_ex^2 >= 0``), so the logarithm is safe.
    """
    eta_plus = np.asarray(eta_plus, dtype=float)
    eta_minus = np.asarray(eta_minus, dtype=float)
    out = np.empty_like(eta_plus)

    big = eta_plus > 30.0
    small = ~big

    if np.any(small):
        s = d_plus * np.cosh(eta_plus[small]) - d_minus * np.cos(eta_minus[small])
        out[small] = np.arccosh(np.maximum(s, 1.0))

    if np.any(big):
        e1 = np.exp(-eta_plus[big])
        ratio = d_minus / d_plus
        out[big] = (
            eta_plus[big]
            + np.log(d_plus)
            + np.log1p(e1 * e1 - 2.0 * ratio * np.cos(eta_minus[big]) * e1)
        )

    return np.maximum(out, 0.0)


def r2eff_carver_richards(
    nu_cpmg: np.ndarray,
    p_b: float,
    k_ex: float,
    dw: float,
    r2_a: float = 0.0,
    r2_b: Optional[float] = None,
) -> np.ndarray:
    """Exact two-site CPMG dispersion profile (Carver-Richards).

    Valid in all exchange regimes; reduces to :func:`r2eff_luz_meiboom` when
    ``k_ex >> dw`` (asserted in the test suite).
    """
    nu = np.asarray(nu_cpmg, dtype=float)
    r2_b = r2_a if r2_b is None else r2_b
    p_a = 1.0 - p_b

    if k_ex <= 0 or dw == 0.0 or p_b <= 0.0 or p_b >= 1.0:
        return np.full_like(nu, 0.5 * (r2_a + r2_b))

    xi = (r2_a - r2_b) - p_a * k_ex + p_b * k_ex
    psi = xi * xi - dw * dw + 4.0 * p_a * p_b * k_ex * k_ex
    zeta = 2.0 * dw * xi
    root = np.sqrt(psi * psi + zeta * zeta)

    d_plus = 0.5 * (1.0 + (psi + 2.0 * dw * dw) / root)
    d_minus = 0.5 * (-1.0 + (psi + 2.0 * dw * dw) / root)

    tau_cp = 1.0 / (4.0 * np.maximum(nu, 1e-30))
    eta_plus = np.sqrt(2.0) * tau_cp * np.sqrt(max(root + psi, 0.0))
    eta_minus = np.sqrt(2.0) * tau_cp * np.sqrt(max(root - psi, 0.0))

    acosh = _arccosh_cr(d_plus, d_minus, eta_plus, eta_minus)
    return 0.5 * (r2_a + r2_b + k_ex - 2.0 * nu * acosh)


def r2eff_from_fit(
    fit: ExchangeFit,
    nu_cpmg: np.ndarray,
    r2_0: float = 0.0,
    exact: bool = True,
) -> np.ndarray:
    """Synthetic dispersion profile for a fitted residue."""
    if exact:
        return r2_0 + r2eff_carver_richards(nu_cpmg, fit.p_minor, fit.k_ex, fit.dw)
    return r2eff_luz_meiboom(nu_cpmg, fit.p_minor, fit.k_ex, fit.dw, r2_0=r2_0)


def exchange_contrast(fit: ExchangeFit) -> float:
    """``p_A p_B * delta_obs^2`` -- the statistic to rank cryptic sites on.

    Not ``R_ex``, and not the CPMG amplitude.  Both of those divide by ``k_ex``,
    and ``k_ex`` is the one fitted quantity whose value in an MD ensemble cannot
    be transferred to the real system: a 100 ns trajectory only resolves
    exchange above ~1e8 s^-1, while the openings we care about run at 1e3-1e6
    s^-1, so any score containing ``1/k_ex`` ranks on a number the simulation
    got wrong by orders of magnitude.

    What *is* transferable is the contrast: how distinct the alternate state is
    (``delta_obs^2``) and how much of it there is (``p_A p_B``).  The squared
    term dominates, which is what separates a genuine opening from a flexible
    loop that visits two nearly identical states.
    """
    return float(fit.p_minor * (1.0 - fit.p_minor) * fit.delta_obs ** 2)


def nmr_window(k_ex: float) -> str:
    """Which solution-NMR experiment could observe exchange at this rate.

    The forward CPMG model's real job is not scoring -- it is telling an
    experimentalist which measurement to run on a predicted site, and warning
    when no measurement can reach it.
    """
    if k_ex <= 0:
        return "none (no exchange fitted)"
    if k_ex < 1e1:
        return "HDX-MS / real-time NMR"
    if k_ex < 1e3:
        return "CEST / DEST"
    if k_ex < 1e5:
        return "CPMG relaxation dispersion"
    if k_ex < 1e7:
        return "R1rho off-resonance"
    return "none (motionally narrowed; MD-timescale only)"


def dispersion_amplitude(fit: ExchangeFit, nu_lo: float = 25.0,
                         nu_hi: float = 1000.0) -> float:
    """``R2eff(nu_lo) - R2eff(nu_hi)`` -- what a real CPMG experiment would see.

    Use this to decide whether a *predicted* excited state is experimentally
    reachable, not to rank sites: it is identically zero for the fast exchange
    that MD actually samples.  :func:`exchange_contrast` is the ranking
    statistic.
    """
    lo, hi = r2eff_from_fit(fit, np.array([nu_lo, nu_hi]))
    return float(lo - hi)

# --------------------------------------------------------------------------
# What the trajectory is capable of resolving.
# --------------------------------------------------------------------------

def detectability_limit(
    total_time_ns: float,
    p_b: float,
    min_transitions: int = 8,
) -> float:
    """Slowest ``k_ex`` (s^-1) a trajectory of this length can resolve.

    A two-state fit needs the excited state to be entered and left several
    times.  The expected number of transitions in a trajectory of total length
    ``T`` is ``2 T p_A p_B k_ex``, so requiring ``min_transitions`` of them sets
    a floor on the exchange rate:

        k_ex_min = min_transitions / (2 T p_A p_B)

    This is the number that decides whether a *negative* dispersion result
    means anything.  Pooling all three ATLAS replicas of one protein gives
    T = 300 ns, and a 5% excited state then needs k_ex >= ~3e8 s^-1 to be seen
    at all -- whereas cryptic-pocket opening is commonly 1e3-1e6 s^-1.  The
    dispersion channel is therefore expected to be silent for most genuine
    cryptic sites in ATLAS-length data, which is precisely why the tail and
    breathing channels exist rather than being redundant confirmations.
    """
    if not 0.0 < p_b < 1.0 or total_time_ns <= 0:
        return float("inf")
    return float(min_transitions / (2.0 * total_time_ns * 1e-9 * p_b * (1.0 - p_b)))


def is_detectable(fit: "ExchangeFit", total_time_ns: float,
                  min_transitions: int = 8) -> bool:
    """Whether ``fit``'s exchange is fast enough for this trajectory to resolve."""
    if fit.p_minor <= 0.0 or fit.p_minor >= 1.0:
        return False
    return fit.k_ex >= detectability_limit(total_time_ns, fit.p_minor, min_transitions)

