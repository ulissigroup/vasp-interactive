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
ediffg = -1


def test_steps():
    from ase.build import molecule

    """Test if VaspInteractive correctly write inputs
    """
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tempdir:
        # Large ediffg but with auto correction
        with VaspInteractive(
            xc="pbe", ediff=ediff, ediffg=ediffg, directory=tempdir
        ) as calc:
            h2.calc = calc
            dyn = BFGS(h2)
            dyn.run(fmax=fmax)
            # EDIFFG == -1 won't work
            assert dyn.nsteps >= 5

        # Try giving ediffg but don't allow overwrite
        with pytest.raises(ValueError):
            calc = VaspInteractive(
                xc="pbe",
                ediff=ediff,
                ediffg=ediffg,
                directory=tempdir,
                allow_default_param_overwrite=False,
            )
