"""Grid-based buried-cavity detection: the measurement the pipeline was missing.

Every observable this project had -- SASA, RMSF, amide exchange, contact number
-- measures how *exposed* a residue is.  A cryptic pocket is a *cavity*, and a
buried cavity has zero solvent-accessible surface area.  So the primary
observable was blind by construction to the thing being detected, which is why
the whole pipeline scored 0.512 against a 0.562 mean-SASA baseline: it was
ranking surface loops.

This implements LIGSITE's protein-solvent-protein (PSP) criterion (Hendlich et
al. 1997; Huang & Schroeder 2006).  A grid point is *free* if no heavy atom's
van der Waals sphere plus a probe radius reaches it, and it is *enclosed along
an axis* when protein lies on both sides of it along that axis.  Free points
enclosed on two axes line a surface pocket; on all three, a buried cavity.

The enclosure test is the whole trick and it vectorises: cumulative-OR of the
occupancy grid from each end of an axis says whether protein lies before and
after every point at once, so all three axes cost three passes over a boolean
array rather than a ray cast per point.

Cavity volume is then assigned to whichever residues line it, which is the
quantity a cryptic-site label actually refers to.
"""

from __future__ import annotations

import numpy as np

# Heavy-atom van der Waals radii, angstroms.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
       "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
DEFAULT_VDW = 1.70


def _occupancy(xyz: np.ndarray, radii: np.ndarray, origin: np.ndarray,
               shape: tuple, spacing: float, probe: float) -> np.ndarray:
    """Boolean grid: True where an atom's vdW sphere plus probe reaches."""
    occ = np.zeros(shape, dtype=bool)
    rmax = float(radii.max() + probe)
    pad = int(np.ceil(rmax / spacing))
    idx = np.rint((xyz - origin) / spacing).astype(int)
    # Offsets within the largest sphere, tested per atom against its own radius.
    rng = np.arange(-pad, pad + 1)
    dx, dy, dz = np.meshgrid(rng, rng, rng, indexing="ij")
    off = np.stack([dx.ravel(), dy.ravel(), dz.ravel()], axis=1)
    d2 = (off ** 2).sum(1) * spacing ** 2
    for a in range(xyz.shape[0]):
        keep = d2 <= (radii[a] + probe) ** 2
        p = idx[a] + off[keep]
        ok = ((p >= 0) & (p < np.array(shape))).all(1)
        p = p[ok]
        occ[p[:, 0], p[:, 1], p[:, 2]] = True
    return occ


def _enclosure(occ: np.ndarray, span: int = 0) -> np.ndarray:
    """Per-point count of axes with protein on both sides (0-3).

    ``span`` bounds how far along the axis the protein may be, in grid cells.
    With an unbounded scan almost every free point near a globular protein
    counts as enclosed on all three axes, because the infinite ray eventually
    strikes the far side of the molecule: on a 112-residue test case that marked
    98 of 112 residues as lining a cavity, which is not a measurement. A bounded
    span asks the local question instead -- is there protein close by on both
    sides -- which is what distinguishes a pocket from open solvent.
    """
    from scipy.ndimage import maximum_filter1d
    psp = np.zeros(occ.shape, dtype=np.uint8)
    for axis in (0, 1, 2):
        if span > 0:
            w = 2 * span + 1
            fwd = maximum_filter1d(occ, size=w, axis=axis, origin=span, mode="constant")
            rev = maximum_filter1d(occ, size=w, axis=axis, origin=-span, mode="constant")
        else:
            fwd = np.maximum.accumulate(occ, axis=axis)
            rev = np.flip(np.maximum.accumulate(np.flip(occ, axis=axis), axis=axis),
                          axis=axis)
        psp += (fwd & rev).astype(np.uint8)
    return psp


def cavity_grid(xyz: np.ndarray, elements: list, spacing: float = 1.0,
                probe: float = 1.4, min_enclosure: int = 3, span_A: float = 8.0):
    """Cavity grid points for one frame.

    ``min_enclosure`` 3 keeps only points buried on all three axes -- interior
    voids, which is what a cryptic pocket is before it opens.  2 admits surface
    grooves as well and is the right setting for comparing against a pocket
    detector rather than a cavity detector.
    """
    radii = np.array([VDW.get(e.upper(), DEFAULT_VDW) for e in elements])
    lo = xyz.min(0) - (radii.max() + probe + 2.0)
    hi = xyz.max(0) + (radii.max() + probe + 2.0)
    shape = tuple(np.ceil((hi - lo) / spacing).astype(int) + 1)
    occ = _occupancy(xyz, radii, lo, shape, spacing, probe)
    free = ~occ
    psp = _enclosure(occ, span=int(round(span_A / spacing)))
    pocket = free & (psp >= min_enclosure)
    pts = np.argwhere(pocket) * spacing + lo
    return pts, float(spacing ** 3)


def residue_cavity_volume(xyz: np.ndarray, elements: list, res_of_atom: np.ndarray,
                          n_res: int, spacing: float = 1.0, probe: float = 1.4,
                          min_enclosure: int = 3, lining: float = 5.0,
                          span_A: float = 8.0) -> np.ndarray:
    """Cavity volume (A^3) lining each residue in one frame.

    A cavity point is credited to every residue with a heavy atom within
    ``lining`` angstroms, because a pocket is lined by all of them and a
    cryptic-site label names all of them too.
    """
    from scipy.spatial import cKDTree
    pts, cell = cavity_grid(xyz, elements, spacing, probe, min_enclosure, span_A)
    out = np.zeros(n_res)
    if pts.size == 0:
        return out
    tree = cKDTree(xyz)
    for j, near in enumerate(tree.query_ball_point(pts, lining)):
        if near:
            for r in set(res_of_atom[near]):
                out[r] += cell
    return out
