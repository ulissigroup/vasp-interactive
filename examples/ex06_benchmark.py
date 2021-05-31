"""Run benchmark for VaspInteractive on different molecule systems
Benchmark molecules taken from https://wiki.fysik.dtu.dk/gpaw/devel/ase_optimize/ase_optimize.html
"""
import numpy as np
import os
import tempfile
import random
from time import time
from pathlib import Path
import tempfile

from ase.db import connect
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from vasp_interactive.vasp_interactive import parse_outcar_iterations
from ase.calculators.vasp import Vasp

curdir = Path(__file__).parent 
systems = []
with connect(curdir / "systems.db") as conn:
    for row in conn.select():
        atoms = row.toatoms()
        # make pbc=True only for VASP
        atoms.pbc = True
        systems.append(atoms)

# default parameters that shared by all vasp calculators
default_params = dict(xc="pbe", ismear=0, sigma=0.01, kspacing=0.5, kgamma=True, npar=4)
fmax = 0.05

# Following functions do the relaxation and returns ionic / electronic steps with wall time

def relax_vasp_interactive(atoms):
    """Vasp Interactive"""
    atoms = atoms.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(directory=tmpdir)
        params.update(default_params)
        calc = VaspInteractive(**params)
        t_ = time()
        with calc:
            atoms.calc = calc
            # Suppress output
            dyn = BFGS(atoms)
            dyn.run(fmax=fmax)
        n_ion, n_elec = calc.read_all_iterations()
        t_wall = time() - t_
        e = atoms.get_potential_energy()
    return e, n_ion - 1, n_elec[:-1], t_wall


def relax_vasp_bfgs(atoms):
    """Classic vasp + bfgs, no wave function reloading"""
    atoms = atoms.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(istart=0, ibrion=-1, nsw=0, directory=tmpdir)
        params.update(default_params)
        calc = Vasp(**params)
        atoms.calc = calc
        dyn = BFGS(atoms, logfile=None)
        n_elec = []
        n_ion = 1
        
        # Use manual force threshold in order to read the iterations
        t_ = time()
        f = np.abs(atoms.get_forces()).max()
        n_elec.append(calc.read_number_of_iterations())
        while f > fmax:
            dyn.step()
            n_ion += 1
            f = np.abs(atoms.get_forces()).max()
            n_elec.append(calc.read_number_of_iterations())
        n_elec = np.array(n_elec)
        t_wall = time() - t_
        e = atoms.get_potential_energy()
    return e, n_ion, n_elec, t_wall

def relax_vasp(atoms):
    """Classic vasp"""
    atoms = atoms.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(istart=0, ediffg=-fmax, ibrion=2, nsw=500, directory=tmpdir)
        params.update(default_params)
        calc = Vasp(**params)
        atoms.calc = calc
        t_ = time()
        atoms.get_potential_energy()
        n_ion, n_elec = parse_outcar_iterations(calc.load_file("OUTCAR"))
        t_wall = time() - t_
        e = atoms.get_potential_energy()
    return e, n_ion, n_elec, t_wall


def main():
    import pickle
    res_file = curdir / "benchmark.pkl"
    if res_file.is_file():
        with open(res_file, "rb") as fd:
            results = pickle.load(fd)
    else:
        results = dict()
    for i in range(len(systems)):
        atoms = systems[i]
        name = atoms.get_chemical_formula()
        if name in results.keys():
            print(f"Results for {name} loaded from pickle")
            continue
        res = dict()
        print(f"Relaxation for {name}")
        
        print("\tVasp Interactive...")
        res["vasp-inter"] = relax_vasp_interactive(atoms)
        
        print("\tVasp BFGS...")
        res["vasp-bfgs"] = relax_vasp_bfgs(atoms)
        
        print("\tVasp only...")
        res["vasp"] = relax_vasp(atoms)
        
        results[name] = res
        
        # Save at each epoch
        with open(res_file, "wb") as fd:
            pickle.dump(results, fd, protocol=3)
    
    

if __name__ == "__main__":
    main()
