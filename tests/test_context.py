"""Test context methods
"""
import pytest
import numpy as np
import os
import tempfile

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

d = 0.9575
h2_root = Atoms(
    "H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True
)


def test_no_context():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),
            directory=tmpdir,
        )

        h2.calc = calc
        dyn = BFGS(h2)
        # simulate 2 steps
        dyn.step()
        dyn.step()
        # The vasp process should still be running
        assert calc.process.poll() is None
        assert calc.final is False
        
        # Manual finalizing
        calc.finalize()
        assert calc.final is True
        assert calc.process is None
        assert h2.get_potential_energy() == pytest.approx(-6.44217, 1e-4)
    return

def test_context():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # Best practice of VaspInteractive is to use it as ContextManager
        with calc:
            h2.calc = calc
            dyn = BFGS(h2)
            # Simulate steps
            dyn.step()
            dyn.step()
            assert calc.steps == 2
            assert calc.final is False
            assert calc.process.poll() is None
        assert calc.final is True
        assert calc.process is None
        assert h2.get_potential_energy() == pytest.approx(-6.44217, 1e-4)
    return

def test_no_context_with_exception():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # Best practice of VaspInteractive is to use it as ContextManager
        # context manager will exit and exception is handled outside, to calc.__exit__ is handled
        try:
            h2.calc = calc
            dyn = BFGS(h2)
            # Now ASE-BFGS controls the relaxation, not VASP
            # Simulate a single step
            dyn.step()
            pytest.raises(RuntimeError)
        except RuntimeError:
            assert calc.process is None
    return


def test_context_with_exception():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # Best practice of VaspInteractive is to use it as ContextManager
        # context manager will exit and exception is handled outside, to calc.__exit__ is handled
        try:
            with calc:
                h2.calc = calc
                dyn = BFGS(h2)
                # Now ASE-BFGS controls the relaxation, not VASP
                # Simulate a single step
                dyn.step()
                pytest.raises(RuntimeError)
        except RuntimeError:
            assert calc.process is None
    return