"""Running NEB calculations with calculator pausing. This is the same example as ex12_k8s_neb.py but 
the individual calculators are sharing the same CPU resources and only 1 calculator is running at the same time.

Example adapted from https://wiki.fysik.dtu.dk/ase/tutorials/neb/diffusion.html

Sample output:

16 MPI processes, VASP 6.3.0

Sequential normal VASP: Elasp time: 1091.16 s
Sequential VaspInteractive (shared calculator): Elasp time: 2756.08 s
Paused VaspInteractive (individual calculator): Elasp time: 437.10 s
"""
import time
import tempfile
import os
from ase.io import read
from ase.neb import NEB
from ase.optimize import BFGS
from copy import deepcopy, copy
from pathlib import Path
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
import numpy as np

curdir = Path(__file__).parent
initial = read(curdir / "neb-trajs" / "initial.traj")
final = read(curdir / "neb-trajs" / "final.traj")
# Same vasp theoretical level for the initial / final structures
# vasp_params = dict(xc="pbe", kspacing=1.0, encut=150, sigma=0.1, istart=0, lwave=False)
vasp_params = dict(xc="pbe", kspacing=0.5, encut=350, sigma=0.1, istart=0, lwave=False)
# Very rough
# fmax = 1.0
fmax = 0.03


class PausedVaspInteractive(VaspInteractive):
    def calculate(
        self,
        atoms=None,
        properties=["energy"],
        system_changes=["positions", "numbers", "cell"],
    ):
        # _ensure_mpi context ensures only within this block the mpi process is running
        with self._ensure_mpi():
            super().calculate(
                atoms=atoms, properties=properties, system_changes=system_changes
            )


def make_images(n_images=5):
    """Interpolate n_images - 2"""
    images = [initial]
    for i in range(n_images - 2):
        images.append(initial.copy())
    images.append(final)
    return images


def run_neb(images, traj_name, **kwargs):
    neb = NEB(images, **kwargs)
    neb.interpolate()
    t_ = time.perf_counter()
    dyn = BFGS(neb, trajectory=traj_name)
    dyn.run(fmax=fmax)
    elasp = time.perf_counter() - t_
    return elasp


def sequential_vasp(n_images=5):
    print("*" * 50)
    print("Running NEB with sequential VASP code")
    images = make_images(n_images)
    with tempfile.TemporaryDirectory() as tmpdir:
        base_calc = Vasp(directory=tmpdir, **vasp_params)
        for img in images[1 : n_images - 1]:
            img.calc = deepcopy(base_calc)
        elasp = run_neb(images, "neb-vasp-seq.traj")
    energies = np.array([img.get_potential_energy() for img in images])
    print("Energies: ", energies)
    print(f"Elasp time: {elasp:.2f} s")
    return images, energies


def sequential_vasp_inter(n_images=5):
    print("*" * 50)
    print("Running NEB with sequential VaspInteractive code")
    images = make_images(n_images)
    with tempfile.TemporaryDirectory() as tmpdir:
        base_calc = VaspInteractive(directory=tmpdir, **vasp_params)
        for img in images[1 : n_images - 1]:
            img.calc = base_calc
        # All calculations share the same calculator instance
        elasp = run_neb(images, "neb-vasp-inter-seq.traj", allow_shared_calculator=True)
        energies = np.array([img.get_potential_energy() for img in images])
        print("Energies: ", energies)
        print(f"Elasp time: {elasp:.2f} s")
        # Finalize
        for img in images[1 : n_images - 1]:
            img.calc.finalize()
    return images, energies


def paused_vasp_inter(n_images=5):
    print("*" * 50)
    print("Running NEB with Pausable VaspInteractive code")
    images = make_images(n_images)
    with tempfile.TemporaryDirectory() as tmpdir:
        # base_calc = VaspInteractive(directory=tmpdir, **vasp_params)
        tmpdir = Path(tmpdir)
        for i, img in enumerate(images[1 : n_images - 1]):
            calc = PausedVaspInteractive(
                directory=tmpdir / f"img{i:03d}", **vasp_params
            )
            img.calc = calc

        # All calculations share the same calculator instance
        elasp = run_neb(images, "neb-vasp-inter-paused.traj")
        energies = np.array([img.get_potential_energy() for img in images])
        print("Energies: ", energies)
        print(f"Elasp time: {elasp:.2f} s")
        # Finalize
        for img in images[1 : n_images - 1]:
            img.calc.finalize()
    return images, energies


if __name__ == "__main__":
    """Please run the code by
    ```
    python ex14_neb_pause.py 2>/dev/null
    ```
    to plunge stderr message
    """
    traj1, eng1 = sequential_vasp(5)
    traj2, eng2 = sequential_vasp_inter(5)
    traj3, eng3 = paused_vasp_inter(5)
    print(eng2 - eng1)
    print(eng3 - eng1)
    # k8s_vasp_inter(5)
