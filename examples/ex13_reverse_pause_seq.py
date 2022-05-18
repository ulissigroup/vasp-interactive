"""Reverse calculator pausing sequence as compared with ex13_pause_mpi.py
In the main thread, the mpi process is paused until explicitly calling the calc.calculate() method
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


def make_single_calc_call(atoms, i):
    """Only run single calculator call using atoms.calc"""
    atoms.calc._resume_calc()
    atoms.rattle(stdev=0.01, seed=random_seeds[i])
    t_vpi_start = time.time()
    e = atoms.get_potential_energy()
    t_vpi_end = time.time()
    print(f"DFT: loop {i}, {t_vpi_end - t_vpi_start} s.")
    atoms.calc._pause_calc()
    return


def main():
    print(f"Testing VaspInteractive with reversed pausing sequence")
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **vasp_params)
        atoms = mol.copy()
        atoms.calc = calc
        atoms.calc._pause_calc()
        print("Calculator paused")
        with calc:
            for i, seed in enumerate(random_seeds[:2]):
                print("Running expensive function")
                multiproc_expensive_function()
                print("Running DFT")
                make_single_calc_call(atoms, i)
    return


if __name__ == "__main__":
    main()
