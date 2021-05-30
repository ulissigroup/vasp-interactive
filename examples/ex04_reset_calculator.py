"""Example showing how to apply VaspInteractive calculators for different atoms
"""
import numpy as np
import os
import tempfile

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

# Methane and ethylene molecules
_m1 = molecule("CH4", vacuum=5, pbc=True)
_m2 = molecule("C2H4", vacuum=5, pbc=True)


def _dyn_run(atoms, calc, fmax=0.05):
    """run dynamics on atoms"""
    atoms.calc = calc
    dyn = BFGS(atoms)
    dyn.run(fmax=fmax)
    return


def using_context():
    """Using VaspInteractive in context mode will automatically
    reset the calculator
    """
    print("*" * 40)
    print("Running VaspInteractive in context mode")
    print("*" * 40)

    m1 = _m1.copy()
    m2 = _m2.copy()

    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),
            directory=tmpdir,
        )
        # Methane
        with calc:
            _dyn_run(m1, calc)
            print(f"Relaxation of methane: {calc.steps} steps")
            print("Final energy: ", m1.get_potential_energy())

        # Ethylene
        with calc:
            _dyn_run(m2, calc)
            print(f"Relaxation of ethylene: {calc.steps} steps")
            print("Final energy: ", m2.get_potential_energy())

    print("")


def no_context():
    """Using VaspInteractive in backward-compatible approach
    In this case user need to finalize and reset the calculator
    """

    m1 = _m1.copy()
    m2 = _m2.copy()

    print("*" * 40)
    print("Running VaspInteractive in backward-compatible mode")
    print("*" * 40)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),
            directory=tmpdir,
        )
        # Methane
        _dyn_run(m1, calc)
        calc.finalize()
        print(f"Relaxation of methane: {calc.steps} steps")
        print("Final energy: ", m1.get_potential_energy())

        # Ethylene
        calc.reset()
        _dyn_run(m2, calc)
        calc.finalize()
        print(f"Relaxation of ethylene: {calc.steps} steps")
        print("Final energy: ", m2.get_potential_energy())

    print("")


if __name__ == "__main__":
    using_context()
    no_context()
