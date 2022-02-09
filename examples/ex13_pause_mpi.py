"""Pausing and resuming VaspInteractive by sending SIGTSTP / SIGCONT signals to mpirun process
This can be expecially useful when there are computationally expensive codes running in between 2
VASP ionic steps. 

Currently such functionality has only been tested in OpenMPI >= 1.3.1. The user may need to enable 
`orte_forward_job_control=1` in the mca setting to make signal-forward working.
"""
import time
import timeit
import tempfile
import os
from ase.io import read
from ase.build import molecule
from copy import deepcopy, copy
from pathlib import Path
import numpy as np
import random
from multiprocessing import Pool
from vasp_interactive import VaspInteractive

curdir = Path(__file__).parent
random_seeds = [133, 357, 274, 331, 140]
vasp_params = dict(xc="pbe", encut=250, istart=0, lwave=False)

mol = molecule("CH3CH2NH2", vacuum=5, pbc=True)


def expensive_function(seed=42):
    """Useless function only to consume cpu"""
    # t_start = time.time()
    size = 1024
    A = np.random.random((size, size))
    B = np.random.random((size, size))
    C = np.dot(A, B)
    # t_end = time.time()
    # t_elasp = t_end - t_start
    t_elasp = timeit.timeit(lambda: np.dot(A, B), number=4)
    return t_elasp


def multiproc_expensive_function(nprocs=8, seed=42):
    with Pool(nprocs) as pool:
        t_list = pool.map(
            expensive_function,
            [
                seed,
            ]
            * nprocs,
        )
    return t_list


def run_calc(pause=False):
    t_start = time.time()
    print(f"Testing VaspInteractive with pause={pause}")
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **vasp_params)
        atoms = mol.copy()
        atoms.calc = calc
        with calc:
            for i, seed in enumerate(random_seeds):
                atoms.rattle(stdev=0.01, seed=seed)
                t_vpi_start = time.time()
                e = atoms.get_potential_energy()
                t_vpi_end = time.time()
                print(f"DFT: loop {i}, {t_vpi_end - t_vpi_start} s.")
                if pause:
                    calc.pause_calc()
                    t_list = multiproc_expensive_function(8, seed)
                    calc.resume_calc()
                else:
                    t_list = multiproc_expensive_function(8, seed)
                print(f"Expensive computation: loop {i}, {np.max(t_list)} s.")
    t_end = time.time()
    print(f"Total computation time {t_end - t_start} s.")
    return


if __name__ == "__main__":
    run_calc(pause=False)
    time.sleep(3)
    run_calc(pause=True)
