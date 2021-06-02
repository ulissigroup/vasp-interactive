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

    # Collect data, only the first 4 systems
    selections = ["H2", "Cu16", "Cu2", "CAu8O"]
    for i in range(len(selections)):
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

    selections = ["H2", "Cu16", "Cu2", "CAu8O"]
    fig, axes = plt.subplots(4, 2, figsize=(10, 12))
    for i, name in enumerate(selections):
        ax1 = axes[i][0]
        ax2 = axes[i][1]

        # N electronic steps
        n1s = []
        n2s = []

        # time
        t1s = []
        t2s = []
        for method in ["GPMin", "BFGS", "QuasiNewton", "FIRE"]:
            e, n, n1, t1 = results[name][method]["vasp-inter"]
            e, n, n2, t2 = results[name][method]["vasp-ase"]
            # number of steps
            n1s.append(np.sum(n1))
            n2s.append(np.sum(n2))
            # time
            t1s.append(t1)
            t2s.append(t2)
        e, n, nrm, trm = results[name]["RMM-DIIS"]
        e, n, ncg, tcg = results[name]["CG"]
        nrm = np.sum(nrm)
        ncg = np.sum(ncg)

        d = np.arange(4)
        w = 0.2
        ax1.bar(d - w, t1s, w * 2, label="VaspInteractive + ASE")
        ax1.bar(d + w, t2s, w * 2, label="Vasp + ASE")
        ax1.bar([4, 5], [trm, tcg], w * 2, label="Pure VASP")
        arr = np.hstack([t1s, t2s, [trm, tcg]])
        ax1.axhline(y=np.min(arr), ls="--", color="grey")
        ax1.set_xticks(np.arange(6))
        ax1.set_xticklabels(["GPMin", "BFGS", "QuasiNewton", "FIRE", "RMM-DIIS", "CG"])
        #         ax1.set_title("")
        ax1.set_ylabel("Wall Time (s)")
        ax1.text(
            x=0.05,
            y=0.95,
            s=f"{name}-Time",
            ha="left",
            va="top",
            transform=ax1.transAxes,
            fontweight="bold",
        )
        #         ax1.legend()

        # steps plot
        ax2.bar(d - w, n1s, w * 2, label="VaspInteractive + ASE")
        ax2.bar(d + w, n2s, w * 2, label="Vasp + ASE")
        ax2.bar([4, 5], [nrm, ncg], w * 2, label="Pure VASP")
        ax2.set_xticks(np.arange(6))
        ax2.set_xticklabels(["GPMin", "BFGS", "QuasiNewton", "FIRE", "RMM-DIIS", "CG"])
        ax2.text(
            x=0.95,
            y=0.95,
            s=f"{name}-SCF",
            ha="right",
            va="top",
            transform=ax2.transAxes,
            fontweight="bold",
        )
        arr = np.hstack([n1s, n2s, [nrm, ncg]])
        ax2.axhline(y=np.min(arr), ls="--", color="grey")
        ax2.set_ylabel(r"Total $N^{\mathrm{SCF}}$")
        ax2.legend(loc=0)

    fig.tight_layout(pad=1.5)
    fig.savefig(curdir / "benchmark-large.png")


# def plot_details(results):
#     import matplotlib.pyplot as plt
#     fig, ax = plt.subplots(1, 1, figsize=(6, 4))
#     name = "CAu8O"
#     disp_name = {"vasp-inter": "VaspInteractive + BFGS",
#                  "vasp-bfgs": "Vasp + BFGS",
#                  "vasp": "Pure VASP"}
#     for met in ("vasp-inter", "vasp-bfgs", "vasp"):
#         steps = results[name][met][2]
#         ax.plot(steps, "-", label=disp_name[met])
#     ax.legend()
#     ax.set_xlabel("Ionic steps")
#     ax.set_ylabel("Electronic SCF per Ionic Cycle")
#     ax.set_title("CO on Au(111) surface (CAu8O)")

#     fig.tight_layout()
#     fig.savefig(curdir / "details.png")

if __name__ == "__main__":
    #     print(globals())
    results = compute()
    plot_benchmark(results)
#     plot_details(results)
