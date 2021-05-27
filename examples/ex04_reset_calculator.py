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
m1 = molecule("CH4", vacuum=5, pbc=True)
m2 = molecule("C2H4", vacuum=5, pbc=True)


def _dyn_run(atoms, calc, fmax=0.05):
    """run dynamics on atoms
    """
    atoms.calc = calc
    dyn = BFGS(atoms)
    dyn.run(fmax=fmax)
    return

def using_context():
    """Using VaspInteractive in context mode will automatically 
    reset the calculator
    """
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
            
        # Ethylene
        with calc:
            print(calc.steps, calc.final)
            _dyn_run(m2, calc)
            print(f"Relaxation of ethylene: {calc.steps} steps")
    

if __name__ == "__main__":
    using_context()
