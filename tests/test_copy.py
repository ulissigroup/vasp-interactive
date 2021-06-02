"""Makesure copy and deepcopy on the calculate work
"""
import pytest
import numpy as np
from vasp_interactive import VaspInteractive
import tempfile
from pathlib import Path
import os
from ase.atoms import Atoms
from ase.optimize import BFGS
from copy import copy, deepcopy

d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
params = dict(xc="pbe", kpts=(1, 1, 1), ismear=0)
rootdir = Path(__file__).parents[1] / "sandbox"
fmax = 0.05
ediff = 1e-4


def test_copy():
    """Test simple copy method"""
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            h2.calc = calc
            new_calc = copy(calc)
            # Before calculation, calc does not have atoms
            assert new_calc.results == {}
            assert new_calc.atoms is None
            del new_calc
            # Perform a dynamic step
            dyn = BFGS(h2)
            dyn.step()
            new_calc = copy(calc)
            # Are the results same as old?
            assert new_calc.results
            # Behavior of shallow copy does not actually copy the inner components
            assert id(new_calc.results) == id(calc.results)
            assert new_calc.get_potential_energy() == calc.get_potential_energy()
            assert np.any(new_calc.get_forces() == calc.get_forces())
            assert np.any(new_calc.get_stress() == calc.get_stress())
            # Only the process is copied
            assert calc.process
            assert new_calc.process is None
            old_e = new_calc.get_potential_energy()
            # Do a second step, shallow copy will leak the calc.results dict
            dyn.step()
            assert new_calc.get_potential_energy() != old_e
            assert new_calc.get_potential_energy() == calc.get_potential_energy()


def test_deepcopy():
    """Test simple copy method"""
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            h2.calc = calc
            new_calc = deepcopy(calc)
            # Before calculation, calc does not have atoms
            assert new_calc.results == {}
            assert new_calc.atoms is None
            del new_calc
            # Perform a dynamic step
            dyn = BFGS(h2)
            dyn.step()
            new_calc = deepcopy(calc)
            # Are the results same as old?
            assert id(new_calc.results) != id(calc.results)
            assert new_calc.get_potential_energy() == calc.get_potential_energy()
            assert np.any(new_calc.get_forces() == calc.get_forces())
            assert np.any(new_calc.get_stress() == calc.get_stress())
            assert calc.process
            assert new_calc.process is None
            old_e = new_calc.get_potential_energy()
            # Do a second step, is the energy maintained in the new_calc?
            dyn.step()
            assert new_calc.get_potential_energy() == old_e
            assert new_calc.get_potential_energy() != calc.get_potential_energy()
