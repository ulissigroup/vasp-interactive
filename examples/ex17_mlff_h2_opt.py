"""Use VaspInteractive to calculate relaxation of H2 molecule using machine learning force field (MLFF)
Following methods are tested:
1. VaspInteractive interfacing Pure DFT
2. Use previously MLFF run dataset to do online MLFF optimization
3. Do a 200-step MD and reuse MLFF for optimization (inference-only)
4. Reuse MLFF from 2000-step MD for optimization (inference-only)

In methods 2-4, VaspInteractive can leverage the built-in MLFF for optimization purposes other than MD.
"""
import numpy as np
import os
import shutil
import tempfile

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp
import os
from pathlib import Path
from packaging.version import Version

curdir = Path(__file__).parent
h2_root = molecule("H2", pbc=True, cell=[8, 8, 8])
h2_root.rattle(0.1)


def dft():
    h2 = h2_root.copy()
    print("1. Running relaxation of H2 molecule using pure DFT")
    with tempfile.TemporaryDirectory() as tempdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            istart=0,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tempdir,
        )
        with calc:
            h2.calc = calc
            opt = BFGS(h2)
            opt.run(fmax=0.01)
        e = h2.get_potential_energy()
        d = h2.get_distance(0, 1)
        print(calc.version)
        if Version(calc.version) < Version("6.3.0"):
            raise Exception("Must run the example with VASP>= 6.3.0")
    return e, d


def mlff_online_opt():
    h2 = h2_root.copy()
    print(
        "2. Running relaxation of H2 molecule using online MLFF with previous dataset"
    )
    with tempfile.TemporaryDirectory() as tempdir:
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            istart=0,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            custom={"ml_lmlff": True, "ml_istart": 1},
            directory=tempdir,
        )
        td = Path(tempdir)
        ab_file = curdir / "mlff" / "ML_AB.H2"
        shutil.copy(ab_file, td / "ML_AB")
        with calc:
            h2.calc = calc
            opt = BFGS(h2)
            opt.run(fmax=0.01)
        e = h2.get_potential_energy()
        d = h2.get_distance(0, 1)
    return e, d


def mlff_train_inference():
    with tempfile.TemporaryDirectory() as tempdir:
        calc = Vasp(
            ismear=0,
            xc="pbe",
            istart=0,
            ibrion=0,
            mdalgo=2,
            tebeg=500,
            potim=1.0,
            nsw=200,  # Increase this number for better convergence
            kpts=(1, 1, 1),  # not important, just keeps it faster
            custom={"ml_lmlff": True, "ml_istart": 0},
            directory=tempdir,
        )
        h2 = h2_root.copy()
        h2.calc = calc
        print("Generating MLFF. May take a few minutes...")
        h2.get_potential_energy()
        td = Path(tempdir)
        print("Copy ML force field checkpoint")
        h2.get_potential_energy()
        shutil.copy(td / "ML_FFN", td / "ML_FF")
        print("3. Running relaxation with reusable MLFF from 200 step MD")
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            istart=0,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            custom={"ml_lmlff": True, "ml_istart": 2},
            directory=tempdir,
        )

        h2 = h2_root.copy()
        h2.calc = calc
        with calc:
            h2.calc = calc
            opt = BFGS(h2)
            opt.run(fmax=0.01)
        e = h2.get_potential_energy()
        d = h2.get_distance(0, 1)
    return e, d


def mlff_copy_checkpoint():
    with tempfile.TemporaryDirectory() as tempdir:
        td = Path(tempdir)
        ff_file = curdir / "mlff" / "ML_FF.H2"
        print("Copy ML force field checkpoint from 2000 step MD")
        h2 = h2_root.copy()
        h2.get_potential_energy()
        shutil.copy(ff_file, td / "ML_FF")
        print("4. Running relaxation with reusable MLFF from 2000 step MD")
        calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            istart=0,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            custom={"ml_lmlff": True, "ml_istart": 2},
            directory=tempdir,
        )
        h2 = h2_root.copy()
        h2.calc = calc
        with calc:
            h2.calc = calc
            opt = BFGS(h2)
            opt.run(fmax=0.01)
        e = h2.get_potential_energy()
        d = h2.get_distance(0, 1)
    return e, d


if __name__ == "__main__":
    e1, d1 = dft()
    e2, d2 = mlff_online_opt()
    e3, d3 = mlff_train_inference()
    e4, d4 = mlff_train_inference()
    print("\tDFT\tOnline\tInference 1\tInterference 2")
    print("Energy:", e1, e2, e3, e4)
    print("H-H distance:", d1, d2, d3, d4)
