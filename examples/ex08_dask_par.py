"""Modified example of ex04 using Dask parallelism
To run this example you need to have `dask` and `dask-kubernetes` installed
Modify the `worker-cpu-spec.yml` file in the same directory and
adapt `PYTHONPATH`, `volumeMounts` and `volumes` to your situation. 

Ideally, you should see the wall time on the master process to be less than 
the sum of time in individual Dask pods. Overheads include creating kubernetes pods
and inter-pod communication.
"""
import numpy as np
import os
import tempfile
import dask
from pathlib import Path
from dask_kubernetes import KubeCluster
from dask.distributed import Client
import dask.bag as dg
from time import time

from ase.build import molecule
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

# Series of molecules
_m1 = molecule("CH4", vacuum=5, pbc=True)
_m2 = molecule("C2H4", vacuum=5, pbc=True)
_m3 = molecule("C6H6", vacuum=5, pbc=True)
_m4 = molecule("CH3COCH3", vacuum=5, pbc=True)



def _dyn_run(atoms):
    """Run dynamics on atoms in Dask, return the relaxed energy and time"""
    fmax = 0.05
    _t = time()
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms.calc.set(directory=tmpdir)
        dyn = BFGS(atoms)
        dyn.run(fmax=fmax)
        atoms.calc.finalize()
        n_ion_scf, _ = atoms.calc.read_all_iterations()
        t_w = time() - _t
        e = atoms.get_potential_energy()

    return e, t_w, n_ion_scf


def run_all():
    """Using Dask parallelism to split workloads
    """
    # Init dask
    curdir = Path(__file__).parent
    yml_file = (curdir / "worker-cpu-spec.yml").as_posix()
    cluster = KubeCluster(yml_file)
    client = Client(cluster)
    cluster.adapt(minimum=0, maximum=4)
    
    
    calc = VaspInteractive(
            ismear=0,
            xc="pbe",
            kpts=(1, 1, 1),
        )
    seq = [_m1.copy(), _m2.copy(), _m3.copy(), _m4.copy()]
    for atoms in seq:
        atoms.calc = calc

    t_ = time()
    seq_computed = dg.from_sequence(seq).map(_dyn_run).compute()
    t_wall = time() - t_
    
    # Print results
    print("System\tTime (s)\tIonic steps")
    for i in range(len(seq)):
        atoms = seq[i]
        e, t, N = seq_computed[i]
        print(f"{atoms.get_chemical_formula()}\t{t:.2f}\t{N}")
    print(f"System wall time (root process): {t_wall:.2f} s")


if __name__ == "__main__":
    run_all()
