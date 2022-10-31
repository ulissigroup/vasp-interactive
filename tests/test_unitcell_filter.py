"""Test the performance of ASE's unit cell filter on VaspInteractive
"""
import ase
import numpy as np
import tempfile
from ase.constraints import UnitCellFilter, ExpCellFilter, StrainFilter
from ase.optimize import BFGS
from ase.io import read, write
from ase.build import bulk
from ase.calculators.vasp import Vasp

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible

# Expanded lattice than optimal
a_ = 2.95
al = bulk("Al", a=a_, b=a_, c=a_)
# Very bad parameter only for speed
params = {
    "xc": "pbe",
    "encut": 200,
    "ediff": 1.0e-4,
    "kpts": [3, 3, 3],
    "gamma": True,
    "isym": 0,
}


def test_uc_filter():
    skip_lattice_if_incompatible()
    # Normal VASP
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = al.copy()
        atoms.calc = Vasp(directory=tmpdir, ibrion=-1, nsw=1, **params)
        # lwave=True, istart=1, **params)
        sf1 = UnitCellFilter(atoms, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf1)
        dyn.run(fmax=0.05)
    # VaspInteractive
    atoms2 = al.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms2.calc = VaspInteractive(directory=tmpdir, **params)
        sf2 = UnitCellFilter(atoms2, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf2)
        dyn.run(fmax=0.05)
        atoms2.calc.finalize()
    print(atoms.cell.cellpar())
    print(atoms2.cell.cellpar())
    # There is still minor discrepancy between normal VASP and VaspInteractive
    assert abs(atoms.cell.cellpar()[0] - atoms2.cell.cellpar()[0]) < 0.005


def test_ec_filter():
    skip_lattice_if_incompatible()
    # Normal VASP
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = al.copy()
        atoms.calc = Vasp(directory=tmpdir, ibrion=-1, nsw=1, **params)
        # lwave=True, istart=1, **params)
        sf1 = ExpCellFilter(atoms, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf1)
        dyn.run(fmax=0.05)
    # VaspInteractive
    atoms2 = al.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms2.calc = VaspInteractive(directory=tmpdir, **params)
        sf2 = ExpCellFilter(atoms2, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf2)
        dyn.run(fmax=0.05)
        atoms2.calc.finalize()
    # There is still minor discrepancy between normal VASP and VaspInteractive
    print(atoms.cell.cellpar())
    print(atoms2.cell.cellpar())
    assert abs(atoms.cell.cellpar()[0] - atoms2.cell.cellpar()[0]) < 0.005


def test_strain_filter():
    skip_lattice_if_incompatible()
    # Normal VASP
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = al.copy()
        atoms.calc = Vasp(directory=tmpdir, ibrion=-1, nsw=1, **params)
        # lwave=True, istart=1, **params)
        sf1 = StrainFilter(atoms, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf1)
        dyn.run(fmax=0.05)
    # VaspInteractive
    atoms2 = al.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms2.calc = VaspInteractive(directory=tmpdir, **params)
        sf2 = StrainFilter(atoms2, mask=[True, True, True, False, False, False])
        dyn = BFGS(sf2)
        dyn.run(fmax=0.05)
        atoms2.calc.finalize()
    print(atoms.cell.cellpar())
    print(atoms2.cell.cellpar())
    # There is still minor discrepancy between normal VASP and VaspInteractive
    assert abs(atoms.cell.cellpar()[0] - atoms2.cell.cellpar()[0]) < 0.005
