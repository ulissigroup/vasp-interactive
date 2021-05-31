"""Use VaspInteractive to calculate relaxation of H2 molecule
"""
import numpy as np
import os
import shutil
import tempfile

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive


def run_with_vasp_interactive():
    print("Running relaxation of H2 molecule using VaspInteractive")
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )
        d = 0.9575
        h2 = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
        # Best practice of VaspInteractive is to use it as ContextManager
        with calc:
            h2.calc = calc
            dyn = BFGS(h2)
            # Now ASE-BFGS controls the relaxation, not VASP
            dyn.run(fmax=0.05)
        # Read the ionic and electronic steps
        # last ionic iteration is dummy, omit
        n_ion, n_elec = calc.read_all_iterations()
        print(f"Relaxation by VaspInteractive in {n_ion - 1} steps")
        print(f"Electronic scf steps per ionic cycle: {n_elec[:-1]}")


if __name__ == "__main__":
    run_with_vasp_interactive()
