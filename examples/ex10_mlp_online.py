"""Running online al_mlp learner using VaspInteractive
   Dependency:
   - https://github.com/ulissigroup/al_mlp
   - https://github.com/ulissigroup/cluster_mlp
   This example requires CPU cores >= 8.0
"""
from al_mlp.atomistic_methods import replay_trajectory
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from ase.io import Trajectory
import numpy as np
import time
from ase.optimize import BFGS
from pathlib import Path
import torch
import os
import copy
import contextlib
import io

default_flare_config = {
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

default_learner_params = {
    "filename": "relax_example",
    "file_dir": "mlp_examples/",
    "stat_uncertain_tol": 0.08,
    "dyn_uncertain_tol": 0.1,
    "fmax_verify_threshold": 0.05,  # eV/AA
}

vasp_flags = {
    "isif": 0,
    "lreal": "Auto",
    "encut": 350.0,
    "ncore": 8,
    "xc": "PBE",
    "txt": "-",
}

def gen_online_calc(images, parent_calc, 
                    flare_config=default_flare_config, 
                    learner_params=default_learner_params):
    """Use default parameters to generate online calc
    """
    from al_mlp.ml_potentials.flare_pp_calc import FlarePPCalc
    from al_mlp.online_learner.online_learner import OnlineLearner
    ml_potential = FlarePPCalc(flare_config, images)
    calc = OnlineLearner(learner_params, images, ml_potential, parent_calc)
    return calc

def gen_cluster(metal="Cu", number=10):
    """Use random code in cluster_mlp to generate a cluster
    """
    from cluster_mlp.fillPool import fillPool
    from ase.data import atomic_numbers, covalent_radii
    eleNames = [metal]
    eleNums = [number]
    eleRadii = [covalent_radii[atomic_numbers[ele]] for ele in eleNames]
    return fillPool(eleNames, eleNums, eleRadii, None)

# curdir = Path("./").resolve()
# example_dir = curdir / "mlp_benchmark"

def parse_scf(output):
    """Parse scf lines from stdout
    """
    import re
    lines = output.split("\n")
    pat = r"DAV\:\s+([\d]+)"
    prev = -1
    elec_steps = []
    for l in lines:
        m = re.match(pat, l)
        if m:
            nex = int(m[1])
            if nex > prev:
                prev = nex
            else:
                elec_steps.append(prev)
                prev = -1
    elec_steps.append(nex)
    return elec_steps

def run_opt(initial_structure, vasp, optimizer=BFGS, 
            use_al=True, 
            store_wf=True, 
            traj_name="oal_relax.traj"):
    """Choose a backend for vasp or VaspInteractive calculator"""
    os.system("rm -rf WAVECAR INCAR POTCAR POSCAR POTCAR")
    assert vasp.lower() in ("vasp", "vaspinteractive")
    calc_name = vasp.lower()
    images = [initial_structure.copy()]

    parent_calc = VaspInteractive(**vasp_flags)

    if calc_name == "vaspinteractive":
        context = parent_calc
    else:
        context = contextlib.suppress()
        parent_calc.set(ibrion=-1, nsw=0)
        if not store_wf:
            parent_calc.set(istart=0, lwave=False)
    
    t_ = time.time()
    output = io.StringIO()
    # Use context lib to redirect output
    with contextlib.redirect_stdout(output):
        with context:
            if use_al:
                real_calc = gen_online_calc(images, parent_calc)
            else:
                real_calc = parent_calc
            images[0].calc = real_calc
            
            
            dyn = optimizer(images[0], 
                            trajectory=traj_name)
            if use_al:
                dyn.attach(replay_trajectory, 1, images[0].calc, dyn)
            dyn.run(fmax=0.05, steps=1000)
    t_elaps = time.time() - t_
        
    
    output_string = output.getvalue()
    elec_steps = parse_scf(output_string)
    final_image = Trajectory(traj_name)[-1]
    return t_elaps, elec_steps, final_image, output_string


if __name__ == "__main__":
    initial_structure = gen_cluster("Cu", 7)
    print("*" * 40)
    print("Running with BFGS + vasp -- no cache")
    t, steps, fin, _ = run_opt(initial_structure,
                          "VASP", 
                          use_al=False, 
                          store_wf=False, )
    print(f"Time: {t:.4s}")
    
    print("*" * 40)
    print("Running with BFGS + vasp -- cache")
    t, steps, fin, _ = run_opt(initial_structure,
                          "VASP", 
                          use_al=False, 
                          store_wf=True, )
    print(f"Time: {t:.4s}")

    print("*" * 40)
    print("Running with BFGS + vasp inter")
    t, steps, fin, _ = run_opt(initial_structure,
                          "VaspInteractive", 
                          use_al=False, )
    print(f"Time: {t:.4s}")



    print("*" * 40)
    print("Running with OAL + vasp")
    t, steps, fin, _ = run_opt(initial_structure,
                          "VASP", 
                          use_al=True, 
                          store_wf=True, )
    print(f"Time: {t:.4s}")
    

    print("*" * 40)
    print("Running with OAL + vasp inter")
    t, steps, fin, _ = run_opt(initial_structure,
                          "VaspInteractive", 
                          use_al=True, )
    print(f"Time: {t:.4s}")

