"""Compare results for relaxation of H2 molecule
   using VaspInteractive vs Vasp internal vs Vasp SinglePoint + BFGS
   using the same force stop criteria the energy results should be almost identical
"""
import numpy as np
import os
import tempfile

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8])


def run_no_context():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # No context manager
        h2.set_calculator(calc)
        dyn = BFGS(h2)
        dyn.run(fmax=0.05)
        # At this point h2 should have finished running
        # but process still active
        if calc.process.poll() is None:
            print("Process still running: ", h2.process)
        else:
            print("Process finished: ", h2.process)


def run_with_exception():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # Best practice of VaspInteractive is to use it as ContextManager
        try:
            with calc:
                h2.set_calculator(calc)
                dyn = BFGS(h2)
                # Now ASE-BFGS controls the relaxation, not VASP
                dyn.run(fmax=0.05)
                print("Final energy using VaspInteractive: ", h2.get_potential_energy())
                raise RuntimeError("Simulate error")
        except Exception as e:
            print("Encountered error", e)
            if calc.process.poll() is None:
                print("Process still running: ", h2.process)
            else:
                print("Process finished: ", h2.process)


if __name__ == "__main__":
    #     run_no_context()
    run_with_exception()
