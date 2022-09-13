"""Check whether VaspInteractive correctly fallback to singlepoint when nsw=0 or 1
"""
import pytest
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from ase.build import molecule
import tempfile
from pathlib import Path
import os
import sys
import random

params = dict(xc="pbe", kpts=(1, 1, 1), nsw=0, ibrion=-1, ismear=0, ediff=1e-2, encut=120)


def test_nsw():
    atoms = molecule("C2H2", vacuum=5, pbc=True)

    with tempfile.TemporaryDirectory() as tempdir:
        atoms0 = atoms.copy()
        atoms0.calc = Vasp(directory=tempdir, **params)
        e0 = atoms0.get_potential_energy()

    with tempfile.TemporaryDirectory() as tempdir:
        atoms1 = atoms.copy()
        with VaspInteractive(directory=tempdir, **params) as calc1:
            atoms1.calc = calc1
            e1 = atoms1.get_potential_energy()

    with tempfile.TemporaryDirectory() as tempdir:
        atoms2 = atoms.copy()
        params["nsw"] = 1
        with VaspInteractive(directory=tempdir, **params) as calc2:
            atoms2.calc = calc2
            e2 = atoms2.get_potential_energy()

    assert e1 == pytest.approx(e0)
    assert e2 == pytest.approx(e0)
    return


def test_mock_sp():
    """Mocking sp calculator using VaspInteractive on molecular sequences"""
    atoms_root = molecule("C2H2", vacuum=5, pbc=True)
    atoms_seq = [atoms_root.copy() for i in range(3)]
    [a.rattle(stdev=0.1, seed=random.randint(1, 1000)) for a in atoms_seq]

    for atoms in atoms_seq:
        # Normal VASP
        with tempfile.TemporaryDirectory() as tempdir:
            atoms0 = atoms.copy()
            atoms0.calc = Vasp(directory=tempdir, **params)
            e0 = atoms0.get_potential_energy()

        with tempfile.TemporaryDirectory() as tempdir:
            atoms1 = atoms.copy()
            with VaspInteractive(directory=tempdir, **params) as calc1:
                atoms1.calc = calc1
                e1 = atoms1.get_potential_energy()

        assert e1 == pytest.approx(e0)
    return
