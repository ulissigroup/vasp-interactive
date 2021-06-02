"""Additional benchmark on various optimization methods
Methods from ASE: 
 - GPmin
 - BFGS
 - QuasiNewton (BFGSLinSearch)
 - FIRE

Methods provided VASP
 - RMM-DIIS (IBRION=1)
 - Conjugate Gradient (IBRION=2)
 
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
from ase.optimize import BFGS, QuasiNewton, GPMin, FIRE
from ase.optimize.sciopt import SciPyFminCG

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

def select_optimizer(method):
    """Use the method name to get optimizer. `method` is a string e.g. `GPMin`"""
    try:
        opt = globals()[method]
    except KeyError:
        opt = None
    return opt

def patched_step(obj, opt):
    """Monkey patch of the step method to add record for the steps"""
    opt.step(obj)
    if not hasattr(obj, "n_elec_scf"):
        obj.n_elec_scf = []
    n_elec = obj.atoms.calc.read_number_of_iterations()
    obj.n_elec_scf.append(n_elec)
    return

# Following functions do the relaxation and returns ionic / electronic steps with wall time


def relax_vasp_interactive(atoms, method="BFGS"):
    """Vasp Interactive"""
    atoms = atoms.copy()
    opt = select_optimizer(method)
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(directory=tmpdir)
        params.update(default_params)
        calc = VaspInteractive(**params)
        t_ = time()
        with calc:
            atoms.calc = calc
            # Suppress output
            dyn = opt(atoms)
            dyn.run(fmax=fmax)
        n_ion, n_elec = calc.read_all_iterations()
        t_wall = time() - t_
        e = atoms.get_potential_energy()
    return e, n_ion - 1, n_elec[:-1], t_wall


def relax_vasp_ase(atoms, method="BFGS"):
    """Classic vasp + ase, no wave function reloading"""
    atoms = atoms.copy()
    opt = select_optimizer(method)
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(istart=0, ibrion=-1, nsw=0, directory=tmpdir)
        params.update(default_params)
        calc = Vasp(**params)
        atoms.calc = calc
        dyn = opt(atoms)
        dyn.step = lambda: patched_step(dyn, opt)
        # ASE dynamics does not record first ionic step, manually simulate
        atoms.get_potential_energy()
        n_elec = []
        n_elec.append(calc.read_number_of_iterations())
        calc.reset()
        # Do actual relaxation
        t_ = time()
        dyn.run(fmax=fmax)
        n_ion = dyn.nsteps + 1
        n_elec.extend(dyn.n_elec_scf)
        n_elec = np.array(n_elec)
        t_wall = time() - t_
        print(n_ion, n_elec)
        e = atoms.get_potential_energy()
    return e, n_ion, n_elec, t_wall

def relax_vasp(atoms, method):
    assert method in ("RMM-DIIS", "CG")
    """Classic vasp"""
    atoms = atoms.copy()
    if method == "RMM-DIIS":
        ibrion = 1
    else:
        ibrion = 2
    with tempfile.TemporaryDirectory() as tmpdir:
        params = dict(istart=0, ediffg=-fmax, ibrion=ibrion, nsw=500, directory=tmpdir)
        params.update(default_params)
        calc = Vasp(**params)
        atoms.calc = calc
        t_ = time()
        atoms.get_potential_energy()
        n_ion, n_elec = parse_outcar_iterations(calc.load_file("OUTCAR"))
        t_wall = time() - t_
        e = atoms.get_potential_energy()
    return e, n_ion, n_elec, t_wall


def compute():
    import pickle
    res_file = curdir / "benchmark-large.pkl"
    if res_file.is_file():
        with open(res_file, "rb") as fd:
            results = pickle.load(fd)
    else:
        results = dict()

    # Collect data
    for i in range(len(systems)):
        atoms = systems[i]
        name = atoms.get_chemical_formula()
        if name in results.keys():
            print(f"Results for {name} loaded from pickle")
        else:
            res = dict()
            print(f"Relaxation for {name}")
            for method in ("GPMin", "BFGS", "QuasiNewton", "FIRE"):
                print(f"\tUsing method {method}")
                par_res = dict()
               
                print("\t\tVasp Interactive...")
                par_res["vasp-inter"] = relax_vasp_interactive(atoms, method)

                print("\t\tVasp ASE...")
                par_res["vasp-ase"] = relax_vasp_ase(atoms, method)
                res[method] = par_res

            for method in ("RMM-DIIS", "CG"):
                print(f"\tUsing Vasp {method}...")
                res[method] = relax_vasp(atoms, method)

            results[name] = res

            # Save at each epoch
            with open(res_file, "wb") as fd:
                pickle.dump(results, fd, protocol=3)

    return results

def plot_benchmark(results):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax1 = axes[0]
    ax2 = axes[1]
    w = 0.15
    # N electronic steps
    n1s = []
    n2s = []
    
    # time
    t1s = []
    t2s = []
    for i, key in enumerate(results.keys()):
        e, n, n1, t1 = results[key]["vasp-inter"]
        e, n, n2, t2 = results[key]["vasp-bfgs"]
        e, n, n3, t3 = results[key]["vasp"]
        # number of steps
        n1s.append(np.sum(n1) - np.sum(n3))
        n2s.append(np.sum(n2) - np.sum(n3))
        # time
        t1s.append(t1 / t3)
        t2s.append(t2 / t3)
        
        

    d = np.arange(len(results))
    w1 = 0.2
    ax1.bar(d - w, t1s, w * 2, label="VaspInteractive + BFGS")
    ax1.bar(d + w, t2s, w * 2, label="Vasp + BFGS")
    ax1.axhline(y=1, ls="--", color="grey")
    ax1.set_xticks(d)
    ax1.set_xticklabels(list(results.keys()))
    ax1.set_title("Rel. Time to Pure VASP")
    ax1.set_ylabel(r"$t / t_{\mathrm{VASP}}$")
    ax1.legend()
    
    # steps plot
    ax2.bar(d - w, n1s, w * 2, label="VaspInteractive + BFGS")
    ax2.bar(d + w, n2s, w * 2, label="Vasp + BFGS")
    ax2.set_xticks(d)
    ax2.axhline(y=0, ls="--", color="grey")
    ax2.set_xticklabels(list(results.keys()))
    ax2.set_title("Rel. Total Electronic SCFs to Pure VASP")
    ax2.set_ylabel(r"$N^{\mathrm{SCF}} - N^{\mathrm{SCF}}_{\mathrm{VASP}}$")
    ax2.legend()
    
    fig.tight_layout()
    fig.savefig(curdir / "benchmark.png")


def plot_details(results):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    name = "CAu8O"
    disp_name = {"vasp-inter": "VaspInteractive + BFGS",
                 "vasp-bfgs": "Vasp + BFGS",
                 "vasp": "Pure VASP"}
    for met in ("vasp-inter", "vasp-bfgs", "vasp"):
        steps = results[name][met][2]
        ax.plot(steps, "-", label=disp_name[met])
    ax.legend()
    ax.set_xlabel("Ionic steps")
    ax.set_ylabel("Electronic SCF per Ionic Cycle")
    ax.set_title("CO on Au(111) surface (CAu8O)")
    
    fig.tight_layout()
    fig.savefig(curdir / "details.png")

if __name__ == "__main__":
#     print(globals())
    results = compute()
#     plot_benchmark(results)
#     plot_details(results)
