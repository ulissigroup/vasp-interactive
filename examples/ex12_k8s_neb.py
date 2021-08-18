"""Running NEB calculations using kubernetes-isolated VaspInteractive processes
parallelization of neb calculations are achieved via multithreading on the local pod

Running NEB calculations using kubernetes ensures process isolation, 
and no need to allocate large number of cpus at the same time. Setup is also easier than MPI-based parallelization

Example adapted from https://wiki.fysik.dtu.dk/ase/tutorials/neb/diffusion.html
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
from vasp_interactive.kubernetes import KubeVaspInteractive, create_kube_pods

curdir = Path(__file__).parent
initial = read(curdir / "neb-trajs" / "initial.traj")
final = read(curdir / "neb-trajs" / "final.traj")
# Same vasp theoretical level for the initial / final structures
vasp_params = dict(xc="pbe", kspacing=0.5, encut=350, sigma=0.1, istart=0, lwave=False)
# Very rough
fmax = 0.25


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
    energies = [img.get_potential_energy() for img in images]
    print("Energies: ", energies)
    print(f"Elasp time: {elasp:.2f} s")
    return


def sequential_vasp_inter(n_images=5):
    print("*" * 50)
    print("Running NEB with sequential VaspInteractive code")
    images = make_images(n_images)
    with tempfile.TemporaryDirectory() as tmpdir:
        base_calc = VaspInteractive(directory=tmpdir, **vasp_params)
        for img in images[1 : n_images - 1]:
            img.calc = base_calc
    # All calculations share the same calculator instance
    elasp = run_neb(images, "neb-vasp-inter-kube.traj", allow_shared_calculator=True)
    energies = [img.get_potential_energy() for img in images]
    print("Energies: ", energies)
    print(f"Elasp time: {elasp:.2f} s")
    # Finalize
    for img in images[1 : n_images - 1]:
        img.calc.finalize()
    return


def k8s_vasp_inter(n_images=5):
    # lazy cleanup
    os.system("rm -rf /home/jovyan/data/scartch/*")
    print("*" * 50)
    print("Running NEB with KubeVaspInteractive code")
    scale = n_images - 2
    cluster, worker_pods = create_kube_pods(scale=n_images - 2, cpu=8, memory="4Gi")
    images = make_images(n_images)
    for i in range(n_images - 2):
        images[i + 1].calc = KubeVaspInteractive(
            directory=f"/home/jovyan/data/scartch/neb-{i + 1}",
            pod=worker_pods[i],
            **vasp_params,
        )
    # All calculations share the same calculator instance
    elasp = run_neb(images, "neb-vasp-inter-kube.traj", parallel=True)
    energies = [img.get_potential_energy() for img in images]
    print("Energies: ", energies)
    print(f"Elasp time: {elasp:.2f} s")

    # not necessary to finalize processes as closing the cluster kills vasp
    cluster.close()
    return


if __name__ == "__main__":
    """Please run the code by
    ```
    python ex12_k8s_neb.py 2>/dev/null
    ```
    to plunge stderr message
    """
    sequential_vasp(5)
    sequential_vasp_inter(5)
    k8s_vasp_inter(5)
