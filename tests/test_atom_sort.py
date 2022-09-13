import pytest
from vasp_interactive import VaspInteractive
from ase.calculators.calculator import CalculatorSetupError
import tempfile
from pathlib import Path
import os
from ase.io import read


rootdir = Path(__file__).parents[1] / "sandbox"
curdir = Path(__file__).parent
cluster_root = read(curdir / "TIP4P-2.xyz")
cluster_root.pbc = True
# Make a dummy "sequence"
parameters = dict(
    xc="pbe", ediff=1e-4, prec="norm", algo="Fast", lreal="auto", lwave=False
)
fmax = 0.05


def test_steps():
    from ase.build import molecule

    """Test if VaspInteractive correctly write inputs
    """
    cluster = cluster_root.copy()
    calc_dir = rootdir / "cluster"
    calc = VaspInteractive(directory=calc_dir, **parameters)
    with calc:
        cluster.calc = calc
        e1 = cluster.get_potential_energy()
        print(e1)
        cluster2 = cluster.copy()
        cluster2.calc = calc
        e2 = cluster2.get_potential_energy()
        print(e2)
        # assert e1 == e2
        cluster3 = cluster.copy()[calc.sort]
        cluster3.calc = calc
        # print(cluster3)
        # # Cannot switch atoms order!
        with pytest.raises((RuntimeError, NotImplementedError)):
            e3 = cluster3.get_potential_energy()
        #     assert e3 == e2
    return
