"""Synthetic hydrogen-deuterium exchange: the cheap, falsifiable channel.

Relaxation dispersion is the most information-rich NMR readout of a cryptic
pocket, but it is not the only one, and it is not the cheapest to reproduce in
silico.  Amide exchange reports the same underlying physics -- transient local
opening -- through a quantity that a phenomenological model predicts well from
structure alone: the protection factor.

The signal to extract is not the protection factor itself but its *anomaly*.
A residue that looks buried and hydrogen-bonded in the deposited structure, yet
exchanges freely across the ensemble, is a residue whose burial is a
time-average over an opening.  That is the structural definition of a cryptic
site, and unlike a dispersion prediction it can be checked in an afternoon by
HDX-MS on any protein, without isotope labelling or a cryoprobe.

The protection-factor model is Best & Vendruscolo's:

    ln PF_i = beta_c <N_c,i> + beta_h <N_h,i>

with ``N_c`` the heavy-atom contact count around the amide nitrogen and
``N_h`` the number of hydrogen bonds accepted by the amide hydrogen.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BETA_C = 0.35
BETA_H = 2.00


@dataclass
class AmideEnvironment:
    """Per-frame burial and hydrogen bonding of each backbone amide."""

    n_contacts: np.ndarray   # (T, R)
    n_hbonds: np.ndarray     # (T, R)


def amide_environment(
    n_coords: np.ndarray,
    h_coords: np.ndarray,
    heavy_coords: np.ndarray,
    heavy_resid: np.ndarray,
    amide_resid: np.ndarray,
    acceptor_mask: np.ndarray,
    contact_cutoff: float = 6.5,
    hbond_cutoff: float = 2.4,
    exclude_window: int = 2,
    chunk: int = 128,
) -> AmideEnvironment:
    """Count contacts and hydrogen bonds per amide, per frame.

    All coordinate arrays are ``(T, ., 3)`` in angstroms.  ``acceptor_mask``
    selects the oxygen/nitrogen acceptors among ``heavy_coords``.

    Neighbouring residues are excluded from the contact count (``exclude_window``)
    because chain connectivity contributes burial that has nothing to do with
    tertiary packing, and including it flattens the very contrast we are after.
    """
    n_coords = np.asarray(n_coords, dtype=float)
    h_coords = np.asarray(h_coords, dtype=float)
    heavy_coords = np.asarray(heavy_coords, dtype=float)
    heavy_resid = np.asarray(heavy_resid, dtype=int)
    amide_resid = np.asarray(amide_resid, dtype=int)
    acceptor_mask = np.asarray(acceptor_mask, dtype=bool)

    t, r, _ = n_coords.shape
    contacts = np.zeros((t, r), dtype=float)
    hbonds = np.zeros((t, r), dtype=float)
    c2 = contact_cutoff ** 2
    h2 = hbond_cutoff ** 2

    acc_idx = np.flatnonzero(acceptor_mask)
    acc_resid = heavy_resid[acc_idx]

    for i in range(r):
        res = amide_resid[i]
        far = np.abs(heavy_resid - res) > exclude_window
        far_acc = np.abs(acc_resid - res) > exclude_window
        if not far.any():
            continue
        for start in range(0, t, chunk):
            stop = min(start + chunk, t)
            d2 = ((heavy_coords[start:stop][:, far, :]
                   - n_coords[start:stop, i, None, :]) ** 2).sum(-1)
            contacts[start:stop, i] = (d2 < c2).sum(axis=1)

            if far_acc.any():
                sel = acc_idx[far_acc]
                dh2 = ((heavy_coords[start:stop][:, sel, :]
                        - h_coords[start:stop, i, None, :]) ** 2).sum(-1)
                hbonds[start:stop, i] = (dh2 < h2).sum(axis=1)

    return AmideEnvironment(contacts, np.minimum(hbonds, 1.0))


def ln_protection_factor(
    n_contacts: np.ndarray,
    n_hbonds: np.ndarray,
    beta_c: float = BETA_C,
    beta_h: float = BETA_H,
) -> np.ndarray:
    """Best-Vendruscolo ``ln PF`` from burial and hydrogen bonding."""
    return beta_c * np.asarray(n_contacts, float) + beta_h * np.asarray(n_hbonds, float)


def ensemble_ln_pf(env: AmideEnvironment, **kw) -> np.ndarray:
    """Ensemble ``ln PF``, averaging the *environment* before the model.

    Averaging inside the linear model (rather than averaging exchange rates) is
    what Best and Vendruscolo parameterised, so it is what reproduces their
    coefficients.
    """
    return ln_protection_factor(env.n_contacts.mean(axis=0),
                                env.n_hbonds.mean(axis=0), **kw)


def rate_averaged_ln_pf(env: AmideEnvironment, **kw) -> np.ndarray:
    """Alternative ``ln PF`` averaging the exchange *rate* (``1/PF``) instead.

    Exchange is a rate process, so the observable is the mean of ``1/PF`` over
    the ensemble, not ``1/mean(PF)``.  This version is dominated by the open
    frames and is therefore far more sensitive to a rare opening -- the gap
    between the two is itself a cryptic-site signal, which is what
    :func:`breathing_anomaly` uses.
    """
    per_frame = ln_protection_factor(env.n_contacts, env.n_hbonds, **kw)
    return -np.log(np.mean(np.exp(-per_frame), axis=0))


def breathing_anomaly(env: AmideEnvironment, **kw) -> np.ndarray:
    """``ln PF(mean environment) - ln PF(rate averaged)``, per residue.

    Zero for a residue whose burial is static.  Large and positive for one whose
    average burial hides brief, deep excursions to a solvent-exposed,
    hydrogen-bond-free state -- a breathing site.  Because exchange is a rate
    process, this difference is exactly the error you make by reading burial off
    a single structure, which is why it isolates what a crystal structure hides.
    """
    return ensemble_ln_pf(env, **kw) - rate_averaged_ln_pf(env, **kw)


def static_vs_ensemble(
    ln_pf_static: np.ndarray,
    env: AmideEnvironment,
    **kw,
) -> np.ndarray:
    """How much more protected the deposited structure looks than the ensemble.

    Use when a reference (crystal/AlphaFold) structure is available: it answers
    the practical question a medicinal chemist asks, which is whether the
    structure in the PDB is lying about how closed this site is.
    """
    return np.asarray(ln_pf_static, float) - rate_averaged_ln_pf(env, **kw)
