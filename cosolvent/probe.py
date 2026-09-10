"""Benzene cosolvent probe, built from CHARMM36's own aromatic parameters.

The point of cosolvent MD is that unbiased simulation cannot reach a cryptic
opening. Cryptic pockets open at roughly 1e3-1e6 s^-1 and 100 ns of unbiased MD
samples 1e8 and faster, so the event is simply absent -- which is what this
project measured three separate ways: cavity score does not track apo-to-holo
pocket RMSD, cryptic residues show *less* excursive dynamics than average, and
the trajectories that scored best were the ones that began in the bound state.

A hydrophobic probe does not wait for the rare event. Benzene partitions into
apolar surface defects and proto-pockets and holds them open, lowering the
barrier by orders of magnitude rather than hoping to sample over it. This is the
mechanism behind CryptoSite's original pipeline and behind SWISH.

**No new force field parameters are introduced.** benzene is exactly the
phenylalanine ring: CHARMM36 types ``CA`` (aromatic carbon, -0.115 e) and ``HP``
(aromatic hydrogen, +0.115 e), whose bonded and Lennard-Jones terms are already
in ``charmm36.xml``. Generating GAFF parameters would have required
openff-toolkit, which is conda-only and unavailable here; reusing validated
aromatic types is both simpler and less likely to be wrong.
"""

from __future__ import annotations

import numpy as np

# CHARMM36 already contains benzene as residue BENZ with proper CGenFF types
# (CG2R61 / HGR61, -0.115/+0.115 e). An earlier version of this file defined a
# duplicate template from the phenylalanine ring types, which OpenMM correctly
# rejected as "Multiple non-identical matching templates". Using the shipped
# residue means no new parameters are introduced at all.
BENZ_ATOMS = ["CG", "HG", "CD1", "HD1", "CD2", "HD2",
              "CE1", "HE1", "CE2", "HE2", "CZ", "HZ"]
BENZ_RING = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]      # ring order
BENZ_BONDS = [("CG", "CD1"), ("CD1", "CE1"), ("CE1", "CZ"), ("CZ", "CE2"),
              ("CE2", "CD2"), ("CD2", "CG"), ("CG", "HG"), ("CD1", "HD1"),
              ("CE1", "HE1"), ("CZ", "HZ"), ("CE2", "HE2"), ("CD2", "HD2")]


CC, CH = 0.1375, 0.1080          # nanometres, CHARMM aromatic geometry


def benzene_coordinates(centre=(0.0, 0.0, 0.0), rng=None):
    """Planar hexagon with a random orientation, in nanometres."""
    ang = np.arange(6) * np.pi / 3.0
    ring = np.stack([CC * np.cos(ang), CC * np.sin(ang), np.zeros(6)], axis=1)
    hyd = np.stack([(CC + CH) * np.cos(ang), (CC + CH) * np.sin(ang),
                    np.zeros(6)], axis=1)
    # Interleaved to match BENZ_ATOMS: CG HG CD1 HD1 ... where ring order is
    # CG CD1 CE1 CZ CE2 CD2, so carbon k of BENZ_RING sits at angle k.
    order = {"CG": 0, "CD1": 1, "CE1": 2, "CZ": 3, "CE2": 4, "CD2": 5}
    xyz = np.empty((12, 3))
    for i, nm in enumerate(BENZ_ATOMS):
        k = order[nm if nm.startswith("C") else
                  {"HG": "CG", "HD1": "CD1", "HD2": "CD2",
                   "HE1": "CE1", "HE2": "CE2", "HZ": "CZ"}[nm]]
        xyz[i] = ring[k] if nm.startswith("C") else hyd[k]
    if rng is not None:
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        xyz = xyz @ q.T
    return xyz + np.asarray(centre)


def place_probes(topology, positions, n_probes, box, rng, min_dist=0.45):
    """Add ``n_probes`` benzenes at random points clear of the protein.

    ``min_dist`` (nm) is the closest a probe centre may sit to any existing
    heavy atom. Probes are placed before solvation so that water fills around
    them rather than having to be deleted.
    """
    from openmm.app import Element, Topology
    from scipy.spatial import cKDTree

    pos = np.array(positions.value_in_unit(positions.unit)) if hasattr(positions, "unit") \
        else np.array(positions)
    tree = cKDTree(pos)
    lo = pos.min(0) - 0.2
    hi = pos.max(0) + 0.2
    added, tries = [], 0
    while len(added) < n_probes and tries < n_probes * 500:
        tries += 1
        c = rng.uniform(lo, hi)
        if tree.query(c)[0] < min_dist + 0.3:
            continue
        xyz = benzene_coordinates(c, rng)
        if added:
            allpts = np.concatenate(added)
            if cKDTree(allpts).query(xyz)[0].min() < 0.35:
                continue
        added.append(xyz)

    top = Topology()
    top.setPeriodicBoxVectors(box)
    chain_map = {}
    for ch in topology.chains():
        chain_map[ch] = top.addChain(ch.id)
    atom_map = {}
    for res in topology.residues():
        nr = top.addResidue(res.name, chain_map[res.chain], res.id)
        for a in res.atoms():
            atom_map[a] = top.addAtom(a.name, a.element, nr)
    for b in topology.bonds():
        top.addBond(atom_map[b[0]], atom_map[b[1]])

    pc = top.addChain("P")
    C, H = Element.getBySymbol("C"), Element.getBySymbol("H")
    for _ in range(len(added)):
        r = top.addResidue("BENZ", pc)
        made = {nm: top.addAtom(nm, C if nm[0] == "C" else H, r)
                for nm in BENZ_ATOMS}
        for a, b in BENZ_BONDS:
            top.addBond(made[a], made[b])
    new_pos = np.concatenate([pos] + added) if added else pos
    return top, new_pos, len(added)
