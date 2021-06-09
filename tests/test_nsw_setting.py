import pytest
from vasp_interactive import VaspInteractive
import tempfile
from pathlib import Path
import os
from ase.atoms import Atoms
from ase.optimize import BFGS

d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
rootdir = Path(__file__).parents[1] / "sandbox"
fmax = 0.05
ediff = 1e-4


def test_single():
    from ase.build import molecule

    """simulate an error that nsw is zero
    """
    h2 = h2_root.copy()
    workdir = rootdir / "test_nsw"

    calc = VaspInteractive(xc="pbe", ediff=ediff, directory=workdir, nsw=0, ibrion=-1, istart=0)
    with calc:
        h2.calc = calc
        e1 = h2.get_potential_energy()
        assert calc.process.poll() is None

        # Simulate second step
        h2.rattle()
        # this will cause error
        with pytest.raises(RuntimeError):
            e2 = h2.get_potential_energy()
    
    return

def test_low_nsw():
    from ase.build import molecule

    """Simulate runtime error when nsw is too small. Always stops after steps >= nsw + 1
    """
    h2 = h2_root.copy()
    workdir = rootdir / "test_nsw2"
    # Use arbitrary
    calc = VaspInteractive(xc="pbe", ediff=ediff, directory=workdir, nsw=5, ibrion=-1, istart=0)
    with calc:
        h2.calc = calc
        for i in range(10):
            h2.rattle()
            try:
                e = h2.get_potential_energy()
            except RuntimeError:
                break
    assert i == 6
    del calc
    
    # Use arbitrary
    calc = VaspInteractive(xc="pbe", ediff=ediff, directory=workdir, nsw=1, ibrion=-1, istart=0)
    with calc:
        h2.calc = calc
        for i in range(10):
            h2.rattle()
            try:
                e = h2.get_potential_energy()
            except RuntimeError:
                break
    assert i == 2
    
    return