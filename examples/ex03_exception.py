"""Example showing garbage collection of VaspInteractive.
Ideally, if running with the context mode, the VASP process is automatically cleaned up
when exiting (e.g. Exception raised)
"""
import numpy as np
import os
import tempfile

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)


def run_no_context():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        #     tmpdir = "./"
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        h2.calc = calc
        dyn = BFGS(h2)
        dyn.run(fmax=0.05)
        # At this point h2 should have finished running
        # but process still active
        if calc.process.poll() is None:
            print("Process still running: ", calc.process)
        else:
            print("Process finished: ", calc.process)
    # Garbage collection will delete calc which then calls `__del__`
    # In this case, VASP process will return >0 code since tmpdir is deleted
    # but no orphan processes will be left
    return


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
        # context manager will exit and exception is handled outside, to calc.__exit__ is handled
        with calc:
            h2.calc = calc
            dyn = BFGS(h2)
            # Now ASE-BFGS controls the relaxation, not VASP
            # Simulate a single step
            dyn.step()
            #             dyn.run(fmax=0.05)
            raise RuntimeError(
                (
                    "Simulate error, please check if there are still orphan processes.\n"
                    "You should see this error message after VaspInteractive closes."
                )
            )
    return


if __name__ == "__main__":
    run_no_context()
    run_with_exception()
