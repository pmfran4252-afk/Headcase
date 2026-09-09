"""Build the benchmark cohort, and say how much it can actually resolve.

The first scan used the 13 proteins where a CryptoBench apo entry and an ATLAS
trajectory happen to be the *same PDB entry and chain*.  That produced a null
result which was then over-read: at n=13 the scan had roughly 40% power against
the effect it observed, so it could not have distinguished a real effect from
nothing.  A null at 40% power is not evidence of absence.

Three fixes, in order of how much they matter.

**1. Match on UniProt, not on PDB entry.**  The exact-entry requirement is an
accident of which crystal form each project happened to pick.  ATLAS publishes a
UniProt id per chain and a per-residue ``UnP_num`` column, and SIFTS maps every
PDB chain to UniProt, so a CryptoBench label can be carried onto a trajectory of
the same protein from a different crystal form with no sequence alignment. That
takes the cohort from 13 to 83 candidates -- past the 35 needed for 80% power on
the strongest observable measured so far.

**2. Guard the apo state.**  A substitute crystal form may be the *holo*
structure, whose cryptic pocket is already open at frame 0.  That is not a
cryptic-site detection test.  ATLAS's ``contact_ligand`` column does not help
here: it is a boolean describing the deposited entry, it is True for 85% of
ATLAS, and the simulations themselves are ligand-free (0 HETATM).  So the guard
has to be conformational, and this module applies the checkable part -- exclude
any candidate whose ATLAS entry is a known ligand-bound member of an apo/holo
pair -- and flags the rest for a per-site openness check before use.

**3. Report the trivial baselines, always.**  The first scan reported the
pipeline's AUROC with no floor to compare it against.  Mean SASA alone scores
0.562 and RMSF alone 0.527 on the original 13, against the pipeline's 0.512 --
so the pipeline was losing to a one-line heuristic, and nothing in the output
said so.  ``baseline_auroc`` exists so that cannot happen silently again.
"""

from __future__ import annotations

import collections
import csv
import json
import pathlib
from dataclasses import dataclass, asdict
from typing import Optional

import numpy as np

SIFTS_URL = ("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
             "pdb_chain_uniprot.tsv.gz")
ATLAS_LIST_URL = ("https://www.dsimb.inserm.fr/ATLAS/data/download/distributions/"
                  "2024_11_18_ATLAS_pdb.txt")


@dataclass
class Candidate:
    uniprot: str
    label_entry: str        # CryptoBench apo chain carrying the labels
    atlas_entry: str        # ATLAS chain providing the trajectory
    same_entry: bool        # True when this is an exact-match pair
    apo_guard: str          # "ok" | "known_holo" | "unchecked"

    def as_dict(self) -> dict:
        return asdict(self)


def sifts_map(path: pathlib.Path) -> dict:
    """``(pdb, chain) -> {uniprot: (pdb_beg, sp_beg)}`` for residue-number transfer."""
    out: dict = collections.defaultdict(dict)
    with path.open() as f:
        next(f)
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                out[(row["PDB"].lower(), row["CHAIN"])][row["SP_PRIMARY"]] = (
                    int(row["PDB_BEG"]), int(row["SP_BEG"]))
            except ValueError:
                continue
    return out


def build(labels: dict, atlas_chains: list, sifts: dict,
          known_holo: Optional[set] = None) -> list:
    known_holo = known_holo or set()
    atlas_by_up: dict = collections.defaultdict(list)
    for a in atlas_chains:
        p, c = a.split("_")
        for up in sifts.get((p.lower(), c), {}):
            atlas_by_up[up].append(a)

    out: list = []
    for pdb, lab in labels.items():
        for chain in sorted({x.split("_")[0] for x in lab}):
            for up in sifts.get((pdb.lower(), chain), {}):
                for a in atlas_by_up.get(up, []):
                    same = a.lower() == f"{pdb}_{chain}".lower()
                    guard = ("ok" if same else
                             "known_holo" if a.split("_")[0].lower() in known_holo
                             else "unchecked")
                    out.append(Candidate(up, f"{pdb}_{chain}", a, same, guard))
    return out


def labels_to_uniprot(labels: list, chain: str, pdb: str, uniprot: str,
                      sifts: dict) -> set:
    """CryptoBench author residue numbers -> UniProt residue numbers."""
    off = sifts.get((pdb.lower(), chain), {}).get(uniprot)
    if off is None:
        return set()
    pdb_beg, sp_beg = off
    shift = sp_beg - pdb_beg
    return {int(x.split("_")[1]) + shift for x in labels
            if x.split("_")[0] == chain and x.split("_")[1].lstrip("-").isdigit()}


def uniprot_to_index(corresp_path: pathlib.Path) -> dict:
    """ATLAS ``UnP_num`` -> trajectory residue index. ATLAS did the alignment."""
    rows = [l.split("\t") for l in corresp_path.read_text().splitlines() if l.strip()]
    head, rows = rows[0], rows[1:]
    col = head.index("UnP_num")
    out = {}
    for idx, r in enumerate(rows):
        if len(r) > col:
            try:
                out[int(r[col])] = idx
            except ValueError:
                pass
    return out


# --------------------------------------------------------------------------
# Reporting that the first scan lacked
# --------------------------------------------------------------------------

def baseline_auroc(auroc_fn, mask: np.ndarray, mean_sasa: np.ndarray,
                   rmsf: np.ndarray) -> dict:
    """Trivial one-line baselines, reported next to every pipeline score.

    Not optional.  On the original 13 proteins the pipeline scored 0.512 while
    mean SASA scored 0.562 -- the pipeline lost to a single geometric
    descriptor, and the report said nothing about it because no floor was
    computed.  A score with a ceiling and no floor is uninterpretable, which is
    the same criticism this project levels at other people's benchmarks.
    """
    def a(v):
        o = np.argsort(-v)
        return auroc_fn(v[o], mask[o])
    return {"mean_sasa": a(mean_sasa), "neg_mean_sasa": a(-mean_sasa),
            "rmsf": a(rmsf)}


def detectable_effect(n: int, sd: float = 0.15, power: float = 0.80,
                      alpha: float = 0.05) -> float:
    """Smallest mean AUROC above 0.5 this cohort size could detect.

    Report it beside the result.  The first scan concluded "does not beat
    chance" from n=13 at about 40% power, which the data did not support: it
    could only have supported "underpowered to tell".
    """
    from scipy import stats
    for delta in np.arange(0.005, 0.50, 0.005):
        d = delta / sd
        crit = stats.t.ppf(1 - alpha / 2, n - 1)
        if stats.nct.sf(crit, n - 1, d * np.sqrt(n)) >= power:
            return float(0.5 + delta)
    return float("nan")


def main() -> int:
    here = pathlib.Path(__file__).resolve().parent
    labels = json.loads((here / "cb_labels.json").read_text())
    atlas = [l.strip() for l in (here / "atlas_pdb.txt").read_text().splitlines()
             if l.strip()]
    sifts = sifts_map(here / "sifts.tsv")

    known_holo = set()
    lp = pathlib.Path("/Volumes/QIDHD/Loki/benchmarks/cryptobench/gate4-manifest.json")
    if lp.exists():
        for p in json.loads(lp.read_text())["pairs"]:
            known_holo.add(p["comparison"]["id"].lower())

    cands = build(labels, atlas, sifts, known_holo)
    exact = [c for c in cands if c.same_entry]
    usable = [c for c in cands if c.apo_guard != "known_holo"]
    print(f"cohort candidates      : {len(cands)}")
    print(f"  exact PDB match      : {len(exact)}   (the original scan)")
    print(f"  excluded, known holo : {sum(c.apo_guard == 'known_holo' for c in cands)}")
    print(f"  usable               : {len(usable)}")
    print(f"    of which unchecked : {sum(c.apo_guard == 'unchecked' for c in usable)}"
          f"  <- need a per-site openness check before use")
    for n in (len(exact), len(usable)):
        print(f"\n  at n={n:>3}: 80% power to detect a mean AUROC of "
              f"{detectable_effect(n):.3f} or better")
    (here / "cohort.json").write_text(
        json.dumps([c.as_dict() for c in usable], indent=2))
    print(f"\nwrote {here / 'cohort.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
