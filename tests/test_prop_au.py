import pytest
import numpy as np
import os
import tempfile
import random
from time import time

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp
from ase.io import Trajectory
from pathlib import Path

rootdir = Path(__file__).parents[1] / "sandbox"
curdir = Path(__file__).parent 

# propylene/Au
traj = Trajectory(curdir / "prop_au_true_relax.traj")
mol_seq = []
for i in range(2):
    # Rattle the molecule and insert into the sequence
    _mm = traj[i].copy()
    _mm.pbc = True
    mol_seq.append(_mm)


def seq_run(calc):
    """Run with calculator on current sequence"""
    es = []
    fs = []
    for i, atoms in enumerate(mol_seq):
        atoms.calc = calc
        e = atoms.get_potential_energy()
        fmax = np.max(np.abs(atoms.get_forces()))
        iters = calc.read_number_of_iterations()
        es.append(e)
        fs.append(fmax)
        print(f"{i}\t{e:.4f}\t{fmax:.4f}\t{iters}")
    return es, fs


def run_vasp_classic():
    """Run with classic Vasp calculator"""
    print("*" * 40)
    print("Running classic VASP calculator on molecule sequence")
    print("*" * 40)
    dir1 = rootdir / "vasp-classic"
    calc = Vasp(xc="pbe", lreal="auto", directory=dir1, lwave=False)
    t_ = time()
    es, fs = seq_run(calc)
    print(f"Wall time: {time() - t_:.2f} s")
    return es, fs


def run_vasp_interactive():
    """Run with classic Vasp calculator"""
    print("*" * 40)
    print("Running interactive VASP calculator on molecule sequence")
    print("*" * 40)
    dir1 = rootdir / "vasp-inter"
    calc = VaspInteractive(xc="pbe", lreal="auto", directory=dir1, lwave=False)
    t_ = time()
    with calc:
        es, fs = seq_run(calc)
    print(f"Wall time: {time() - t_:.2f} s")
    return es, fs


def test_main():
    ec, fc = run_vasp_classic()
    ei, fi = run_vasp_interactive()
    print(ec, fc)
    print(ei, fi)
    for i in range(2):
        assert ec[i] == pytest.approx(ei[i], 0.01)
        assert fc[i] == pytest.approx(fi[i], 0.01)

    