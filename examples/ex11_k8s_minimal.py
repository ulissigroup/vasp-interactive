"""Simple example showing isolating vasp process in separate kubernetes pod
The example requires to be run within an actual kubernetes pod with `kubectl` access.
Dask-kubernetes serves as the backend for pod scheduling

The example uses shared volume mounting on both local and remote pods, 
so it may not work when setting directory to tmpdir.

Assume the volume is mounted to /home/jovyan/dat
"""
import time
from dask_kubernetes import KubeCluster
from ase.build import molecule
from vasp_interactive.kubernetes import KubeVaspInteractive, create_kube_pods
from pathlib import Path
from threading import Thread


mol1 = molecule("CH4", vacuum=5, pbc=True)
mol2 = mol1.copy()
mol2.rattle(stdev=0.1)
vasp_params = dict(xc="pbe", kpts=1, encut=350, istart=0)

root = Path("/home/jovyan/data/scratch")


def _thread_calculate(atoms, energy):
    """A threaded version of atoms.get_potential_energy. Energy is a one-member list
    ideas taken from https://wiki.fysik.dtu.dk/ase/_modules/ase/neb.html#NEB
    """
    energy[0] = atoms.get_potential_energy()
    return


def main():
    # Scheduling and deployment part
    scale = 2
    cluster, worker_pods = create_kube_pods(scale=scale, cpu=8, memory="2Gi")
    print(f"Created a KubeCluster with scale of {scale}...")
    image = cluster.workers[0]._pod.spec.containers[0].image
    print(f"Pods running on image {image}")
    for number, pod in worker_pods.items():
        print(f"Pod {number}, name: {pod['name']}, namespace: {pod['namespace']}")

    # calculation part
    calc1 = KubeVaspInteractive(
        directory=root / "kube-vpi-test1", pod=worker_pods[0], **vasp_params
    )
    calc2 = KubeVaspInteractive(
        directory=root / "kube-vpi-test2", pod=worker_pods[1], **vasp_params
    )

    # Sequential
    with calc1, calc2:
        mol1.calc = calc1
        mol2.calc = calc2
        t_ = time.perf_counter()
        e1 = mol1.get_potential_energy()
        e2 = mol2.get_potential_energy()
        print("Sequential mode:")
        print(e1, e2)
        print(f"Walltime for 2 sp calculations: {time.perf_counter() - t_}")

    # Pseudo-parallel
    with calc1, calc2:
        mol1.calc = calc1
        mol2.calc = calc2

        # need to use mutable object to store energy
        e1 = [999]
        e2 = [999]
        threads = [
            Thread(target=_thread_calculate, args=(mol1, e1)),
            Thread(target=_thread_calculate, args=(mol2, e2)),
        ]
        t_ = time.perf_counter()
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        print("Threaded mode:")
        print(e1[0], e2[0])
        print(f"Walltime for 2 sp calculations: {time.perf_counter() - t_}")

    cluster.close()


if __name__ == "__main__":
    main()
