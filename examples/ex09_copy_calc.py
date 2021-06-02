"""Similar to ex05, but testing the copy / deepcopy functionalities
To reliably copy the parameters of the calculator, use deepcopy() 
instead of copy() to avoid race conditions.

The __deepcopy__ method of VaspInteractive is designed to reset calc.process to None
after copying. 
"""
import numpy as np
import os
import tempfile
import random
from time import time

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp
from copy import deepcopy

# Methane and ethylene molecules
_m1 = molecule("C2H2", vacuum=5, pbc=True)
mol_seq1 = []
for i in range(6):
    # Rattle the molecule and insert into the sequence
    _mm = _m1.copy()
    _mm.rattle(stdev=0.1, seed=random.randint(1, 1000))
    mol_seq1.append(_mm)

_m2 = molecule("C2H4", vacuum=5, pbc=True)
mol_seq2 = []
for i in range(6):
    # Rattle the molecule and insert into the sequence
    _mm = _m2.copy()
    _mm.rattle(stdev=0.1, seed=random.randint(1, 1000))
    mol_seq2.append(_mm)


def seq_run(seq, calc):
    """Run with calculator on current sequence, copy the VaspInteractive calculator and restart"""
    new_calc = deepcopy(calc)
    with tempfile.TemporaryDirectory() as tmp:
        new_calc.set(directory=tmp)
        with new_calc:
            for i, atoms in enumerate(seq):
                atoms.calc = new_calc
                e = atoms.get_potential_energy()
                iters = new_calc.read_number_of_iterations()
                print(f"{i}\t{e:.4f}\t{iters}")
    return


def run_vasp_interactive():
    """Run with classic Vasp calculator"""
    print("*" * 40)
    print("Running interactive VASP calculator on molecule sequence")
    print("*" * 40)
    base_calc = VaspInteractive(xc="pbe")
    print("Tesing relaxation on C2H2")
    seq_run(mol_seq1, base_calc)
    print("Tesing relaxation on C2H2")
    seq_run(mol_seq2, base_calc)
    return


if __name__ == "__main__":
    run_vasp_interactive()
