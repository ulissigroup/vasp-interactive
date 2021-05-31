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


def test_steps():
    from ase.build import molecule

    """Test if VaspInteractive correctly write inputs
    """
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        calc = VaspInteractive(xc="pbe", ediff=ediff, directory=tempdir)
        h2.calc = calc
        # Initialization
        assert calc.steps == 0
        dyn = BFGS(h2)
        dyn.run(fmax=fmax)
        # Dynamics does not count the 1 step
        assert dyn.nsteps == 6
        assert dyn.nsteps == (calc.steps - 1)

        # Manually finalize the calculator
        calc.finalize()

        # calc.steps does not count the extra dummy ionic scf
        print(calc.steps)

        assert calc.final is True
        assert dyn.nsteps == (calc.steps - 1)
        n_ion_scf, n_elec_scf = calc.read_all_iterations()

        # From OUTCAR there are 8 ionic steps
        assert calc.steps == n_ion_scf - 1

        # Last ionic step must be only 1: a dummy step for saving wavefunction
        assert n_elec_scf[-1] == 1
    return
