"""Running online al_mlp learner using VaspInteractive
   Need https://github.com/ulissigroup/al_mlp as dependency.
   This example requires CPU >= 8.0
"""
from al_mlp.atomistic_methods import Relaxation, replay_trajectory
from al_mlp.online_learner.online_learner import OnlineLearner
from al_mlp.ml_potentials.flare_pp_calc import FlarePPCalc
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from ase.io import Trajectory
import numpy as np
from ase.cluster.decahedron import Decahedron
from ase.optimize import BFGS
import torch
import os
import copy

import contextlib


from pathlib import Path

curdir = Path("./").resolve()
example_dir = curdir / "mlp_examples"

os.environ[
    "VASP_COMMAND"
] = "mpirun -q -np 8 --mca btl_vader_single_copy_mechanism none --mca mpi_cuda_support 0 /opt/vasp.6.1.2_pgi_mkl/bin/vasp_gam"

# Initialize with a Decahedron
initial_structure = Decahedron("Cu", 2, 1, 0)
initial_structure.rattle(0.1)
initial_structure.set_pbc(True)
initial_structure.set_cell([15, 15, 15])
initial_structure.center()
name = initial_structure.symbols

# OAL_example configs
flare_config = {
    "sigma": 4.5,
    "power": 2,
    "cutoff_function": "quadratic",
    "cutoff": 5.0,
    "radial_basis": "chebyshev",
    "cutoff_hyps": [],
    "sigma_e": 0.009,
    "sigma_f": 0.005,
    "sigma_s": 0.0006,
    "hpo_max_iterations": 50,
    "freeze_hyps": 0,
}

learner_params = {
    "filename": "relax_example",
    "file_dir": "mlp_examples/",
    "stat_uncertain_tol": 0.08,
    "dyn_uncertain_tol": 0.1,
    "fmax_verify_threshold": 0.05,  # eV/AA
}

vasp_flags = {
    #         "ibrion": -1,
    #         "nsw": 0,
    "isif": 0,
    "isym": 0,
    "lreal": "Auto",
    "symprec": 1e-10,
    "encut": 350.0,
    "ncore": 8,
    "xc": "PBE",
    "txt": "-",
}


def run_opt(vasp, optimizer=BFGS, use_al=True, store_wf=True, traj_name="oal.traj"):
    """Choose a backend for vasp or VaspInteractive calculator"""
    assert vasp in (Vasp, VaspInteractive)
    calc_name = vasp.name

    images = [initial_structure.copy()]
    elements = np.unique(initial_structure.get_chemical_symbols())

    parent_calc = vasp(**vasp_flags)
    calc_dir = example_dir / f"{calc_name}_inter_tmp"
    parent_calc.set(directory=calc_dir)
    os.system(f"rm -rf {calc_dir.as_posix()}/*")

    if calc_name == "VaspInteractive":
        context = parent_calc
    else:
        # For normal vasp use a dummy context
        context = contextlib.suppress()
        parent_calc.set(ibrion=-1, nsw=0)
        if not store_wf:
            parent_calc.set(istart=0, lwave=False)

    ml_potential = FlarePPCalc(flare_config, images)

    with context:
        if use_al:
            real_calc = OnlineLearner(learner_params, images, ml_potential, parent_calc)
        else:
            real_calc = parent_calc
        images[0].calc = real_calc
        dyn = optimizer(images[0], trajectory=(example_dir / traj_name).as_posix())
        if use_al:
            dyn.attach(replay_trajectory, 1, images[0].calc, dyn)
        dyn.run(fmax=0.05, steps=1000)

    #     print(onlinecalc.parent_calls)
    return


if __name__ == "__main__":
    import time

    times = []
    print("*" * 40)
    print("Running with BFGS + vasp -- no cache")
    # force sync of output
    time.sleep(5)
    t_ = time.time()
    run_opt(Vasp, use_al=False, store_wf=False, traj_name="vasp_bfgs_nocache.traj")
    time.sleep(5)
    t_ = time.time() - t_
    print(f"wall time: {t_} s")
    times.append(t_)

    print("*" * 40)
    print("Running with BFGS + vasp -- cache")
    time.sleep(5)
    t_ = time.time()
    run_opt(Vasp, use_al=False, store_wf=True, traj_name="vasp_bfgs_cache.traj")
    t_ = time.time() - t_
    time.sleep(5)
    print(f"wall time: {t_} s")
    times.append(t_)

    print("*" * 40)
    print("Running with BFGS + vasp inter")
    time.sleep(5)
    t_ = time.time()
    run_opt(VaspInteractive, use_al=False, traj_name="vpi_bfgs.traj")
    t_ = time.time() - t_
    time.sleep(5)
    print(f"wall time: {t_} s")
    times.append(t_)

    print("*" * 40)
    print("Running with OAL + vasp")
    time.sleep(5)
    t_ = time.time()
    run_opt(Vasp, use_al=True, store_wf=True, traj_name="vasp_al.traj")
    t_ = time.time() - t_
    time.sleep(5)
    print(f"wall time: {t_} s")
    times.append(t_)

    print("*" * 40)
    print("Running with OAL + vasp inter")
    time.sleep(5)
    t_ = time.time()
    run_opt(VaspInteractive, use_al=True, store_wf=True, traj_name="vpi_al.traj")
    t_ = time.time() - t_
    time.sleep(5)
    print(f"wall time: {t_} s")
    times.append(t_)

    np.save(example_dir / "times.npy", times)


#     run_opt(Vasp)
#     run_opt(VaspInteractive, use_al=False)
