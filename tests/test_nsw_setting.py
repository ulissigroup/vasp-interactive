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

    """for nsw=0, the calculation can still continue
    """
    h2 = h2_root.copy()
    workdir = rootdir / "test_nsw"

    calc = VaspInteractive(
        xc="pbe", ediff=ediff, directory=workdir, nsw=0, ibrion=-1, istart=0
    )
    with calc:
        h2.calc = calc
        e1 = h2.get_potential_energy()
        # at this stage calculation should have finished
        assert calc.process.poll() == 0
        pid1 = calc.process.pid
        print(calc.process)

        # Simulate second step, a new vasp process will start
        h2.rattle()
        e2 = h2.get_potential_energy()
        pid2 = calc.process.pid
        assert pid1 != pid2
        assert e1 != e2

    return


def test_low_nsw():
    from ase.build import molecule

    """Simulate runtime error when nsw is too small. Always stops after steps >= nsw + 1
    """
    h2 = h2_root.copy()
    workdir = rootdir / "test_nsw2"
    # Use arbitrary
    calc = VaspInteractive(
        xc="pbe",
        ediff=ediff,
        directory=workdir,
        nsw=5,
        ibrion=-1,
        istart=0,
        allow_restart_process=False,
    )
    with calc:
        h2.calc = calc
        for i in range(10):
            h2.rattle()
            try:
                e = h2.get_potential_energy()
            except RuntimeError:
                break
    assert i == 5
    del calc

    # Use arbitrary
    calc = VaspInteractive(
        xc="pbe",
        ediff=ediff,
        directory=workdir,
        nsw=1,
        ibrion=-1,
        istart=0,
        allow_restart_process=False,
    )
    with calc:
        h2.calc = calc
        for i in range(10):
            h2.rattle()
            try:
                e = h2.get_potential_energy()
            except RuntimeError:
                break
    assert i == 1

    return


def test_restart():
    from ase.build import molecule

    """Simulate auto restart of vasp if nsw reached
    """
    h2 = h2_root.copy()
    workdir = rootdir / "test_nsw2"
    # Use arbitrary
    calc = VaspInteractive(
        xc="pbe",
        ediff=ediff,
        directory=workdir,
        nsw=5,
        ibrion=-1,
        istart=0,
        allow_restart_process=True,
    )
    with calc:
        h2.calc = calc
        for i in range(10):
            h2.rattle()
            try:
                e = h2.get_potential_energy()
            except RuntimeError:
                break
    assert i == 9

    return
