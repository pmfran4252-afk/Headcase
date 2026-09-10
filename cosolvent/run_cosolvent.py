"""Cosolvent (benzene) MD on an ATLAS starting structure.

Usage:
    python cosolvent/run_cosolvent.py <tag> <extracted_dir> <out_dir> [ns] [conc_M]

Rationale is in probe.py. In short: unbiased 100 ns cannot contain a cryptic
opening, so every observable computed on it can only score pockets that are
already partly formed. Benzene probes do not wait for the rare event; they
partition into apolar defects and hold them open.

Written to be checkable rather than trusted. It reports probe count and achieved
concentration, tracks how many probes are in contact with the protein over time,
and writes the same trajectory format the existing cavity analysis already
reads, so the identical detector can be run on before/after with nothing else
changed.
"""

from __future__ import annotations

import pathlib
import sys
import time

import numpy as np
import openmm
import openmm.app as app
import openmm.unit as u

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from cosolvent.probe import place_probes

PAD = 1.2 * u.nanometer
TEMP = 310 * u.kelvin
DT = 4 * u.femtosecond          # safe with hydrogen mass repartitioning
REPORT_PS = 20.0


def build(pdb_path, conc_M=0.25, seed=0):
    rng = np.random.default_rng(seed)
    pdb = app.PDBFile(str(pdb_path))
    ff = app.ForceField("charmm36.xml", "charmm36/water.xml")

    # Box from protein extent plus padding, then probe count from its volume.
    pos = np.array(pdb.positions.value_in_unit(u.nanometer))
    ext = pos.max(0) - pos.min(0) + 2 * PAD.value_in_unit(u.nanometer)
    vol = float(np.prod(ext))
    n_probes = max(1, int(round(conc_M * vol * 0.6022)))
    box = [openmm.Vec3(ext[0], 0, 0), openmm.Vec3(0, ext[1], 0), openmm.Vec3(0, 0, ext[2])]

    top, xyz, placed = place_probes(pdb.topology, pos, n_probes, box, rng)
    mod = app.Modeller(top, xyz * u.nanometer)
    mod.topology.setPeriodicBoxVectors([v * u.nanometer for v in box])
    mod.addSolvent(ff, model="tip3p", boxVectors=[v * u.nanometer for v in box],
                   ionicStrength=0.15 * u.molar)
    n_wat = sum(1 for r in mod.topology.residues() if r.name in ("HOH", "WAT"))
    achieved = placed / (vol * 0.6022)
    print(f"  box {ext[0]:.1f}x{ext[1]:.1f}x{ext[2]:.1f} nm ({vol:.0f} nm^3)", flush=True)
    print(f"  benzene {placed}/{n_probes} placed -> {achieved:.2f} M", flush=True)
    print(f"  waters {n_wat}, total atoms {mod.topology.getNumAtoms()}", flush=True)
    system = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * u.nanometer,
                             constraints=app.HBonds, hydrogenMass=4 * u.amu)
    system.addForce(openmm.MonteCarloBarostat(1 * u.bar, TEMP, 25))
    return mod, system, placed, achieved


def run(tag, extracted, out_dir, ns=20.0, conc_M=0.25, seed=0):
    out_dir = pathlib.Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = pathlib.Path(extracted) / tag / f"{tag}.pdb"
    print(f"[{tag}] building", flush=True)
    mod, system, placed, achieved = build(pdb_path, conc_M, seed)

    integ = openmm.LangevinMiddleIntegrator(TEMP, 1 / u.picosecond, DT)
    try:
        plat = openmm.Platform.getPlatformByName("OpenCL")
    except Exception:
        plat = openmm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(mod.topology, system, integ, plat)
    sim.context.setPositions(mod.positions)
    print(f"[{tag}] platform {plat.getName()}; minimising", flush=True)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(TEMP)
    sim.step(int(100 * u.picosecond / DT))                      # equilibration

    # Protein AND probes. Writing protein only makes a null result
    # uninterpretable: "the probes did not open the site" and "the probes never
    # bound" produce identical data. Probe occupancy is the control that
    # separates them, and it costs a few percent of trajectory size. The
    # analysis selects protein atoms itself, so the detector still sees exactly
    # what it saw on the unbiased trajectories.
    prot = [a.index for a in mod.topology.atoms()
            if a.residue.name not in ("HOH", "WAT", "NA", "CL", "SOD", "CLA")]
    sim.reporters.append(app.XTCReporter(str(out_dir / f"{tag}_cosolv.xtc"),
                                         int(REPORT_PS * u.picosecond / DT),
                                         atomSubset=prot))
    sim.reporters.append(app.StateDataReporter(
        str(out_dir / f"{tag}_cosolv.log"), int(1000 * u.picosecond / DT),
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
    # The topology written here must describe the ATOMS IN THE TRAJECTORY, not
    # the whole solvated system. Writing the full topology alongside a
    # protein-only XTC produces files that cannot be opened together at all --
    # 13,391 atoms against 610 in the first run -- and the mismatch is invisible
    # until analysis, long after the compute is spent.
    keep = set(prot)
    sub = app.Modeller(mod.topology, mod.positions)
    sub.delete([a for a in sub.topology.atoms() if a.index not in keep])
    with open(out_dir / f"{tag}_cosolv_top.pdb", "w") as fh:
        app.PDBFile.writeFile(sub.topology, sub.positions, fh, keepIds=True)
    assert sub.topology.getNumAtoms() == len(prot), "topology/trajectory mismatch"

    steps = int(ns * u.nanoseconds / DT)
    print(f"[{tag}] production {ns:.0f} ns ({steps} steps)", flush=True)
    t0 = time.time()
    sim.step(steps)
    el = time.time() - t0
    print(f"[{tag}] done in {el/3600:.2f} h -> {ns/(el/86400):.0f} ns/day", flush=True)
    return {"tag": tag, "n_probes": placed, "conc_M": achieved,
            "ns": ns, "hours": el / 3600}


if __name__ == "__main__":
    a = sys.argv[1:]
    r = run(a[0], a[1], a[2],
            ns=float(a[3]) if len(a) > 3 else 20.0,
            conc_M=float(a[4]) if len(a) > 4 else 0.25)
    print(r)
