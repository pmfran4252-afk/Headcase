"""Fetch and prepare an apo structure from the PDB for cosolvent MD.

Needed because the occupancy cohort is no longer limited to proteins with a
pre-computed ATLAS trajectory. Cosolvent MD makes its own trajectory, so any
deposited apo structure qualifies -- which took the eligible pool from 28
chains to 775 fresh UniProt clusters.

Preparation is deliberately conservative. Only the labelled chain is kept, only
ATOM records and the first altloc, and waters and heteroatoms are dropped.
PDBFixer then rebuilds missing heavy atoms and any gaps short enough to model;
an entry needing more repair than that is rejected rather than silently
patched, because a rebuilt loop is an invention and this pipeline scores
geometry.
"""

from __future__ import annotations

import pathlib
import urllib.request

import openmm.app as app
import openmm.unit as u

MAX_REBUILT_RESIDUES = 8


def fetch(pdb_id: str, cache: pathlib.Path) -> pathlib.Path:
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{pdb_id.lower()}.pdb"
    if not p.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        with urllib.request.urlopen(url, timeout=180) as r:
            p.write_bytes(r.read())
    return p


def isolate_chain(src: pathlib.Path, chain: str, out: pathlib.Path) -> pathlib.Path:
    keep = [l for l in src.read_text().splitlines()
            if l.startswith("ATOM") and l[21] == chain and l[16] in (" ", "A")]
    if not keep:
        raise ValueError(f"chain {chain} has no ATOM records in {src.name}")
    out.write_text("\n".join(keep) + "\nEND\n")
    return out


def prepare(pdb_id: str, chain: str, work: pathlib.Path):
    """Return ``(topology, positions, n_rebuilt)`` ready for solvation."""
    from pdbfixer import PDBFixer
    raw = fetch(pdb_id, work / "raw")
    iso = isolate_chain(raw, chain, work / f"{pdb_id.lower()}_{chain}.pdb")

    fixer = PDBFixer(filename=str(iso))
    fixer.findMissingResidues()
    n_rebuilt = sum(len(v) for v in fixer.missingResidues.values())
    if n_rebuilt > MAX_REBUILT_RESIDUES:
        # Terminal gaps are harmless to drop; interior ones would be invented.
        raise ValueError(f"{pdb_id}_{chain} needs {n_rebuilt} rebuilt residues "
                         f"(limit {MAX_REBUILT_RESIDUES})")
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions, n_rebuilt
