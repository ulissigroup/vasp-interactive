import pytest
from vasp_interactive import VaspInteractive
from ase.calculators.calculator import CalculatorSetupError
import tempfile
from pathlib import Path
import os
from ase.atoms import Atoms

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
        try:
            h2.get_potential_energy()
        except CalculatorSetupError:
            assert calc.version[0] == "5"
    return
