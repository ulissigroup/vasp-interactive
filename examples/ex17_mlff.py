"""Use VaspInteractive to calculate relaxation of H2 molecule
"""
import numpy as np
import os
import shutil
import tempfile

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
import os
from pathlib import Path

curdir = Path(__file__).parent



def main():
    print("Running relaxation of H2 molecule using MLFF inference")
    calc = VaspInteractive(
        ismear=0,
        xc="pbe",
        kpts=(1, 1, 1),  # not important, just keeps it faster
        directory="h2_mlff",
        custom={"ml_lmlff": True, "ml_istart": 2}
    )
    # Create the directory first
    calc._ensure_directory()
    # Copy the ML_FF into directory
    shutil.copy(curdir / "mlff" / "ML_FF.H2", Path(calc.directory) / "ML_FF")
    h2 = molecule("H2", pbc=True, cell=[8, 8, 8])
    h2.calc = calc
    with calc:
        # breakpoint()
        h2.get_potential_energy()
        h2.get_forces()
        h2.get_stress()

        h2.rattle(0.05)
        h2.get_potential_energy()
        h2.get_forces()
        h2.get_stress()
    return


if __name__ == "__main__":
    main()
