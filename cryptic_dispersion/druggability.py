"""Stage 3: is the opened site worth a ligand?

Deliberately separate from the openability score, and never blended into it.
Openability and ligandability are different claims with different failure modes:
a site can open wide and hold nothing, or be beautifully enclosed and never
open. Folding them into one number reproduces the error that closed Gate 5 in
the sister project -- two claims in one endpoint, so neither can be scored
alone -- and it is also, empirically, what degraded every blended score tried
on the tuning cohort.

The other thing this module fixes is *when* the question is asked. An earlier
attempt weighted apo cavity volume by lining hydrophobicity and scored 0.442,
below chance. That is the right question asked of the wrong object: it tested
whether a *closed*, water-sized void is hydrophobic. Druggability is a property
of the open state, so every measurement here is taken on the frames where the
site is actually open.

Three properties decide it, following the fpocket/DoGSite line of work:

``volume``      a fragment needs roughly 100 A^3 and a lead-sized ligand 300;
                below ~150 there is nothing to bind.
``enclosure``   mean protein-solvent-protein count of the cavity points, 0-3.
                A shallow surface dimple binds nothing even when it is large.
``apolar``      fraction of lining residues with positive Kyte-Doolittle
                hydropathy. Buried polar cavities are usually water sites.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np

from .cavity import VDW, DEFAULT_VDW, _enclosure, _occupancy

KD = {"ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
      "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "HSD": -3.2,
      "HSE": -3.2, "HSP": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9,
      "MET": 1.9, "PHE": 2.8, "PRO": -1.6, "SER": -0.8, "THR": -0.7,
      "TRP": -0.9, "TYR": -1.3, "VAL": 4.2}

MIN_FRAGMENT_VOLUME = 150.0   # A^3
MIN_ENCLOSURE = 2.2           # of 3
MIN_APOLAR = 0.35


@dataclass
class Druggability:
    residue: int
    opened_volume: float
    enclosure: float
    apolar_fraction: float
    n_lining: int
    verdict: str
    frame: int = -1

    def as_dict(self) -> dict:
        return asdict(self)


def _components(xyz, elements, spacing, probe, min_enclosure, span_A):
    """Connected cavity components, with their enclosure depth."""
    from scipy import ndimage
    radii = np.array([VDW.get(e.upper(), DEFAULT_VDW) for e in elements])
    lo = xyz.min(0) - (radii.max() + probe + 2.0)
    hi = xyz.max(0) + (radii.max() + probe + 2.0)
    shape = tuple(np.ceil((hi - lo) / spacing).astype(int) + 1)
    occ = _occupancy(xyz, radii, lo, shape, spacing, probe)
    psp = _enclosure(occ, span=int(round(span_A / spacing)))
    pocket = (~occ) & (psp >= min_enclosure)
    lab, n = ndimage.label(pocket)
    return lab, n, psp, lo, spacing


def assess_site(xyz: np.ndarray, elements: list, res_of_atom: np.ndarray,
                res_names: list, residue: int, *, spacing: float = 1.0,
                probe: float = 1.4, min_enclosure: int = 2, span_A: float = 8.0,
                lining: float = 5.0, site_radius: float = 9.0,
                frame: int = -1) -> Druggability:
    """Assess the cavity component adjacent to ``residue`` in one open frame.

    ``min_enclosure`` defaults to 2 here rather than 3: an *opened* cryptic site
    is by definition no longer a sealed interior void, so requiring burial on
    all three axes would reject exactly the state being assessed.

    ``site_radius`` bounds the measurement to cavity within that distance of the
    residue. Connected-component analysis alone does not work here: at probe 1.4
    with enclosure 2 the cavity percolates into one surface-spanning region, and
    the first version of this function reported 1,300-4,000 A^3 pockets with two
    different residues returning an identical 4,060, which is a connected
    surface layer rather than a binding site. A pocket is local, so the site is
    defined locally.
    """
    from scipy.spatial import cKDTree
    lab, n, psp, lo, sp = _components(xyz, elements, spacing, probe,
                                      min_enclosure, span_A)
    if n == 0:
        return Druggability(residue, 0.0, 0.0, 0.0, 0, "no_cavity", frame)

    mine = np.flatnonzero(res_of_atom == residue)
    if mine.size == 0:
        return Druggability(residue, 0.0, 0.0, 0.0, 0, "no_atoms", frame)
    idx = np.argwhere(lab > 0)
    pts = idx * sp + lo
    tree = cKDTree(xyz[mine])
    near = tree.query_ball_point(pts, lining)
    touching = {lab[tuple(i)] for i, nb in zip(idx, near) if nb}
    if not touching:
        return Druggability(residue, 0.0, 0.0, 0.0, 0, "no_cavity", frame)

    # Cavity local to this residue: points in a touching component and within
    # site_radius of one of its atoms.
    keep = np.array([bool(nb) and lab[tuple(i)] in touching
                     for i, nb in zip(idx, tree.query_ball_point(pts, site_radius))])
    if not keep.any():
        return Druggability(residue, 0.0, 0.0, 0.0, 0, "no_cavity", frame)
    sel = idx[keep]
    best_vol = float(keep.sum()) * sp ** 3
    enclosure = float(np.mean([psp[tuple(i)] for i in sel]))

    cav_pts = sel * sp + lo
    all_tree = cKDTree(xyz)
    lining_res = set()
    for nb in all_tree.query_ball_point(cav_pts, lining):
        lining_res.update(res_of_atom[nb].tolist())
    kd = [KD.get(res_names[r], 0.0) for r in sorted(lining_res)]
    apolar = float(np.mean([k > 0 for k in kd])) if kd else 0.0

    if best_vol < MIN_FRAGMENT_VOLUME:
        verdict = "too_small"
    elif enclosure < MIN_ENCLOSURE:
        verdict = "too_open"
    elif apolar < MIN_APOLAR:
        verdict = "too_polar"
    else:
        verdict = "ligandable"
    return Druggability(residue, best_vol, enclosure, apolar, len(lining_res),
                        verdict, frame)


def assess_ranked(traj, cavity_series: np.ndarray, residues, *,
                  top_n: int = 20, **kw) -> list:
    """Assess the top-ranked sites, each on its own most-open frame.

    ``cavity_series`` is the ``(frames, residues)`` cavity volume already
    computed for the openability score, so the open frame costs nothing to find.
    """
    heavy = [a.index for a in traj.topology.atoms if a.element.symbol != "H"]
    el = [traj.topology.atom(i).element.symbol for i in heavy]
    res_of = np.array([traj.topology.atom(i).residue.index for i in heavy])
    names = [r.name for r in traj.topology.residues]
    out = []
    for r in list(residues)[:top_n]:
        f = int(np.argmax(cavity_series[:, r]))
        out.append(assess_site(traj.xyz[f, heavy, :] * 10.0, el, res_of, names,
                               int(r), frame=f, **kw))
    return out
