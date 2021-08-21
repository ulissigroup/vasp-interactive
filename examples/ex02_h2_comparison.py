"""Compare results for relaxation of H2 molecule
   using VaspInteractive vs Vasp internal vs Vasp SinglePoint + BFGS
   using the same force stop criteria the energy results should be almost identical
   
   Even if the WAVECAR is reload during the relaxation, VaspInteractive is still faster 
   than classic Vasp due to shorter spin-up time
"""
import numpy as np
import os
import tempfile
from time import time

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from vasp_interactive.vasp_interactive import parse_outcar_iterations
from ase.calculators.vasp import Vasp
from pathlib import Path

d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
rootdir = Path(__file__).resolve().parents[1] / "sandbox"
fmax = 0.005
ediff = 1e-7


def run_with_vasp_interactive():
    """Relaxation using interactive VASP"""
    print("*" * 40)
    print("Running relaxation using VaspInteractive")
    print("*" * 40)
    h2 = h2_root.copy()
    #     with tempfile.TemporaryDirectory() as tmpdir:
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            istart=0,
            ismear=0,
            ediff=ediff,
            xc="pbe",
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
        )

        # Best practice of VaspInteractive is to use it as ContextManager
        t_ = time()
        with calc:
            h2.calc = calc
            dyn = BFGS(h2)
            # Now ASE-BFGS controls the relaxation, not VASP
            dyn.run(fmax=fmax)
            print(
                "Final energy using VaspInteractive: ",
                h2.get_potential_energy(),
            )
        # Read the iterations, should be outside the context
        # Omit the last ionic step since it's a dummy
        n_ion, n_elec = calc.read_all_iterations()
        print(f"Ionic steps: {n_ion - 1}")
        print(f"Electronic scf per ionic cycle: {n_elec[:-1]}")
        print(f"Wall time: {time() - t_:.2f} s")
        print("")


def run_with_vasp():
    """Relaxation using VASP internal routines"""
    print("*" * 40)
    print("Running relaxation using internal VASP routines")
    print("*" * 40)
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(
            istart=0,
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),
            directory=tmpdir,
            ibrion=2,
            nsw=10,
            ediff=ediff,
            ediffg=-fmax,
        )
        # Best practice of VaspInteractive is to use it as ContextManager
        h2.calc = calc
        t_ = time()
        print(
            "Final energy from Vasp internal routines: ",
            h2.get_potential_energy(),
        )
        # Read the iterations
        n_ion, n_elec = parse_outcar_iterations(calc.load_file("OUTCAR"))
        print(f"Ionic steps: {n_ion}")
        print(f"Electronic scf per ionic cycle: {n_elec}")
        print(f"Wall time: {time() - t_:.2f} s")
        print("")


def run_with_vasp_bfgs():
    """Relaxation using BFGS + classic Vasp calculator"""
    print("*" * 40)
    print("Running relaxation using BFGS + classic Vasp (no WAVECAR reload)")
    print("*" * 40)
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(
            ismear=0,
            istart=0,
            xc="pbe",
            ediff=ediff,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
            ibrion=-1,
            nsw=0,
        )
        t_ = time()
        h2.calc = calc
        dyn = BFGS(h2)
        # Use manual force threshold in order to read the iterations
        n_elec = []
        n_ion = 1
        f = np.abs(h2.get_forces()).max()
        n_elec.append(calc.read_number_of_iterations())
        while f > fmax:
            dyn.step()
            n_ion += 1
            f = np.abs(h2.get_forces()).max()
            n_elec.append(calc.read_number_of_iterations())

        print("Final energy from Vasp + BFGS: ", h2.get_potential_energy())
        # Read the iterations
        print(f"Ionic steps: {n_ion}")
        print(f"Electronic scf per ionic cycle: {n_elec}")
        print(f"Wall time: {time() - t_:.2f} s")
        print("")


def run_with_vasp_bfgs_cache():
    """Relaxation using BFGS + classic Vasp calculator"""
    print("*" * 40)
    print("Running relaxation using BFGS + classic Vasp (with WAVECAR reload)")
    print("*" * 40)
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(
            ismear=0,
            xc="pbe",
            ediff=ediff,
            kpts=(1, 1, 1),  # not important, just keeps it faster
            directory=tmpdir,
            ibrion=-1,
            nsw=0,
        )
        t_ = time()
        h2.calc = calc
        dyn = BFGS(h2)
        # Use manual force threshold in order to read the iterations
        n_elec = []
        n_ion = 1
        f = np.abs(h2.get_forces()).max()
        n_elec.append(calc.read_number_of_iterations())
        while f > fmax:
            dyn.step()
            n_ion += 1
            f = np.abs(h2.get_forces()).max()
            n_elec.append(calc.read_number_of_iterations())

        print("Final energy from Vasp + BFGS: ", h2.get_potential_energy())
        # Read the iterations
        print(f"Ionic steps: {n_ion}")
        print(f"Electronic scf per ionic cycle: {n_elec}")
        print(f"Wall time: {time() - t_:.2f} s")
        print("")


if __name__ == "__main__":
    run_with_vasp_interactive()
    run_with_vasp()
    run_with_vasp_bfgs()
    run_with_vasp_bfgs_cache()
