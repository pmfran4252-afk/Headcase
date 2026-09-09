"""Transfer CryptoBench labels onto an ATLAS trajectory by sequence alignment.

SIFTS cannot do this for most candidates. Its `pdb_chain_uniprot` table leaves
``PDB_BEG`` empty for the majority of entries -- and the PDBe API returns
``author_residue_number: None`` for the same rows -- so there is no recorded
author-to-UniProt offset for the labelled chain. That gap, not trajectory
availability, is what limits the cohort: of the fresh candidates with enough
labels, 21 of 25 clusters were unusable for this reason alone.

Both structures are available, so the mapping does not need UniProt as an
intermediary. The label chain and the ATLAS chain are the same protein, usually
the same sequence, so aligning them directly gives the residue correspondence
and is checkable against the structures themselves.

difflib is sufficient and is used deliberately: these are two crystal forms of
one protein, so the alignment is a near-identity match with indels at the
termini and disordered loops, which longest-matching-block recursion handles
exactly. A gap-penalty aligner would add dependencies and tunable parameters to
a problem that has neither.
"""

from __future__ import annotations

import difflib
import pathlib
import urllib.request

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
}


def fetch_pdb(pdb_id: str, cache: pathlib.Path) -> pathlib.Path:
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{pdb_id.lower()}.pdb"
    if not p.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        with urllib.request.urlopen(url, timeout=120) as r:
            p.write_bytes(r.read())
    return p


def chain_sequence(pdb_path: pathlib.Path, chain: str):
    """``(sequence, [author_resnum])`` for one chain, from ATOM/HETATM records."""
    seq, nums, seen = [], [], set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != chain or line[12:16].strip() != "CA":
            continue
        res = line[17:20].strip()
        if res not in THREE_TO_ONE:
            continue
        key = line[22:27]                       # resSeq + insertion code
        if key in seen:
            continue
        seen.add(key)
        try:
            nums.append(int(line[22:26]))
        except ValueError:
            continue
        seq.append(THREE_TO_ONE[res])
    return "".join(seq), nums


def atlas_sequence(corresp: pathlib.Path):
    """``(sequence, [residue_index])`` for the trajectory, from ATLAS's own table."""
    rows = [l.split("\t") for l in corresp.read_text().splitlines() if l.strip()]
    head, rows = rows[0], rows[1:]
    col = head.index("PDB_seq")
    seq, idx = [], []
    for i, r in enumerate(rows):
        if len(r) > col and r[col].strip():
            seq.append(r[col].strip())
            idx.append(i)
    return "".join(seq), idx


def transfer(label_pdb: str, label_chain: str, corresp: pathlib.Path,
             label_residues: set, cache: pathlib.Path):
    """Map author residue numbers on the label chain to trajectory indices.

    Returns ``(indices, identity, coverage)``. ``identity`` is the aligned
    fraction and should be near 1 for two forms of one protein; a low value
    means the pairing is wrong and the entry must be dropped rather than scored.
    """
    lseq, lnums = chain_sequence(fetch_pdb(label_pdb, cache), label_chain)
    aseq, aidx = atlas_sequence(corresp)
    if not lseq or not aseq:
        return [], 0.0, 0.0
    sm = difflib.SequenceMatcher(None, lseq, aseq, autojunk=False)
    pair = {}
    matched = 0
    for a, b, n in sm.get_matching_blocks():
        matched += n
        for k in range(n):
            pair[lnums[a + k]] = aidx[b + k]
    identity = matched / min(len(lseq), len(aseq))
    hits = [pair[r] for r in sorted(label_residues) if r in pair]
    coverage = len(hits) / max(len(label_residues), 1)
    return hits, identity, coverage
