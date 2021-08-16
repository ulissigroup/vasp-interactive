"""Running online al_mlp learner using VaspInteractive
   Need https://github.com/ulissigroup/al_mlp as dependency.
   This example requires CPU >= 8.0
"""
from al_mlp.utils import compute_with_calc
from al_mlp.atomistic_methods import Relaxation
from al_mlp.online_learner.online_learner import OnlineLearner
from al_mlp.ml_potentials.amptorch_ensemble_calc import AmptorchEnsembleCalc
from amptorch.trainer import AtomsTrainer
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


def run_opt(vasp):
    """Choose a backend for vasp or VaspInteractive calculator"""
    assert vasp in (Vasp, VaspInteractive)
    calc_name = vasp.name

    images = []
    elements = np.unique(initial_structure.get_chemical_symbols())
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
    }

    parent_calc = vasp(**vasp_flags)
    calc_dir = example_dir / f"{calc_name}_inter_tmp"
    parent_calc.set(directory=calc_dir)
    os.system(f"rm -rf {calc_dir.as_posix()}")

    if calc_name == "VaspInteractive":
        context = parent_calc
    else:
        # For normal vasp use a dummy context
        context = contextlib.suppress()
        parent_calc.set(ibrion=-1, nsw=0, istart=0)

    with context:
        # Run relaxation with active learning
        OAL_initial_structure = compute_with_calc(
            [initial_structure.copy()], parent_calc
        )[0]

    print(OAL_initial_structure)

    OAL_relaxation = Relaxation(
        OAL_initial_structure, BFGS, fmax=0.02, steps=200, maxstep=0.04
    )

    Gs = {
        "default": {
            "G2": {
                "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
                "rs_s": [0],
            },
            "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
            "cutoff": 6,
        },
    }

    learner_params = {
        "max_iterations": 10,
        "samples_to_retrain": 1,
        "filename": "relax_example",
        "file_dir": "./",
        "stat_uncertain_tol": 0.1,
        "dyn_uncertain_tol": 0.5,
        "fmax_verify_threshold": 0.05,  # eV/AA
        "relative_variance": True,
        "n_ensembles": 2,
        "use_dask": False,
    }

    config = {
        "model": {"get_forces": True, "num_layers": 3, "num_nodes": 5},
        "optim": {
            "device": "cpu",
            "force_coefficient": 4.0,
            "lr": 1,
            "batch_size": 10,
            "epochs": 100,
            "optimizer": torch.optim.LBFGS,
            "optimizer_args": {"optimizer__line_search_fn": "strong_wolfe"},
        },
        "dataset": {
            "raw_data": images,
            "val_split": 0,
            "elements": elements,
            "fp_params": Gs,
            "save_fps": False,
            "scaling": {"type": "standardize"},
        },
        "cmd": {
            "debug": False,
            "run_dir": "./",
            "seed": 1,
            "identifier": "test",
            "verbose": False,
            # "logger": True,
            "single-threaded": True,
        },
    }

    dbname = (
        "reg_" + str(initial_structure.get_chemical_formula()) + "_oal_" + calc_name
    )
    dbname = (example_dir / dbname).as_posix()
    trainer = AtomsTrainer(config)

    #     print("Parent calc process is: ", parent_calc.process)
    ml_potential = AmptorchEnsembleCalc(trainer, learner_params["n_ensembles"])

    with context:
        online_calc = OnlineLearner(
            learner_params,
            images,
            ml_potential,
            parent_calc,
        )

        real_calc = online_calc

        OAL_relaxation.run(real_calc, filename=dbname)

    OAL_image = OAL_relaxation.get_trajectory(dbname)[-1]

    print(
        "Final Image Results:"
        + "\nEnergy:\n"
        + str(OAL_image.get_potential_energy())
        + "\nForces:\n"
        + str(OAL_image.get_forces())
    )
    print("Steps: ", online_calc.parent_calls, online_calc.parent_electronic_steps)
    np.save(
        example_dir / f"elec_steps_{calc_name}.npy", online_calc.parent_electronic_steps
    )
    return


if __name__ == "__main__":
    run_opt(Vasp)
    run_opt(VaspInteractive)
