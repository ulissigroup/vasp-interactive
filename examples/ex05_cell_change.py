"""Check how vasp interactive can be applied to a list of randomly rattled atoms
   Even if the positions for new atoms are totally random, VaspInteractive still 
   overperforms classic Vasp calculator, due to internal caching of wavefunction.
"""
from pathlib import Path
import numpy as np
import os
import tempfile
import random
from time import time

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp
rootdir = Path(__file__).resolve().parents[1] / "sandbox"

_m = molecule("CH4", vacuum=5, pbc=True)
mol_seq = []
for i in range(6):
    # Rattle the molecule and insert into the sequence
    fac = 1 + (i + 1) / 10
    _mm = _m.copy()
    _mm.rattle(stdev=0.1, seed=random.randint(1, 1000))
    _mm.cell = _mm.cell * [fac, fac, fac]
    print(_mm)
    mol_seq.append(_mm)


def seq_run(calc):
    """Run with calculator on current sequence"""
    for i, atoms in enumerate(mol_seq):
        atoms.calc = calc
        e = atoms.get_potential_energy()
        iters = calc.read_number_of_iterations()
        print(f"{i}\t{e:.4f}\t{iters}")
    return


# def run_vasp_classic():
#     """Run with classic Vasp calculator"""
#     print("*" * 40)
#     print("Running classic VASP calculator on molecule sequence")
#     print("*" * 40)
#     with tempfile.TemporaryDirectory() as tmpdir:
#         calc = Vasp(istart=0, xc="pbe", directory=tmpdir)
#         t_ = time()
#         seq_run(calc)
#         print(f"Wall time: {time() - t_:.2f} s")
#     return


def run_vasp_interactive():
    """Run with classic Vasp calculator"""
    print("*" * 40)
    print("Running interactive VASP calculator on molecule sequence")
    print("*" * 40)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(istart=0, xc="pbe", directory=rootdir / "cell")
        t_ = time()
        with calc:
            seq_run(calc)
        print(f"Wall time: {time() - t_:.2f} s")
    return


if __name__ == "__main__":
    # run_vasp_classic()
    run_vasp_interactive()
