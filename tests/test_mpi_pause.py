"""Testing the cpu percent and pause / resume functionalities
Note the test needs to run VASP with MPI cores > 4. Otherwise it will skip.
"""
import pytest
import numpy as np
from vasp_interactive import VaspInteractive
import psutil
import tempfile
from pathlib import Path
import os
from ase.atoms import Atoms
from ase.optimize import BFGS
from copy import copy, deepcopy
from _common_utils import skip_probe, skip_slurm, get_average_cpu


d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
params = dict(xc="pbe", kpts=(1, 1, 1), ismear=0)
rootdir = Path(__file__).parents[1] / "sandbox"
fmax = 0.05
ediff = 1e-4


# def get_average_cpu(pid, interval=0.5):
#     proc = psutil.Process(pid)
#     vasp_procs = [p for p in proc.children(recursive=True) if "vasp" in p.name()]
#     cpu_per = [p.cpu_percent(interval) for p in vasp_procs]
#     return np.mean(cpu_per)


def test_pause_cpu_percent():
    """Send pause signal to mpi process and see if drops below threshold"""
    skip_probe(4, skip_oversubscribe=True)
    # skip_slurm()
    h2 = h2_root.copy()
    threshold_high = 75.0
    threshold_low = 25.0
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            h2.calc = calc
            h2.get_potential_energy()
            pid = h2.calc.process.pid
            cpu_nonstop = get_average_cpu()
            print(cpu_nonstop)
            assert cpu_nonstop > threshold_high
            h2.calc._pause_calc()
            cpu_stop = get_average_cpu()
            print(cpu_stop)
            assert cpu_stop < threshold_low
            h2.calc._resume_calc()
    return


def test_pause_context():
    """Context mode"""
    skip_probe(4, skip_oversubscribe=True)
    # skip_slurm()
    h2 = h2_root.copy()
    threshold_high = 75.0
    threshold_low = 25.0
    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(directory=tmpdir, **params) as calc:
            h2.calc = calc
            h2.get_potential_energy()
            pid = h2.calc.process.pid
            # Context
            with h2.calc.pause():
                cpu_stop = get_average_cpu()
                print(cpu_stop)
                assert cpu_stop < threshold_low

            cpu_nonstop = get_average_cpu()
            print(cpu_nonstop)
            assert cpu_nonstop > threshold_high
    return


def test_always_true():
    """Placeholder"""
    pass
