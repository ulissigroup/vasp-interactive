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

import signal
from contextlib import contextmanager
from _common_utils import skip_probe, skip_slurm


class TimeoutException(Exception):
    """Simple class for timeout"""

    pass


@contextmanager
def time_limit(seconds):
    """Usage:
    try:
        with time_limit(60):
            do_something()
    except TimeoutException:
        raise
    """

    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)



d = 0.9575
h2_root = Atoms("H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], pbc=True)
params = dict(xc="pbe", kpts=(1, 1, 1), ismear=0)
rootdir = Path(__file__).parents[1] / "sandbox"
fmax = 0.05
ediff = 1e-4


def get_average_cpu(pid, interval=0.5):
    proc = psutil.Process(pid)
    vasp_procs = [p for p in proc.children(recursive=True) if "vasp" in p.name()]
    cpu_per = [p.cpu_percent(interval) for p in vasp_procs]
    return np.mean(cpu_per)


def test_paused_close():
    skip_probe(4)
    skip_slurm()
    """Context mode"""
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        h2.calc = calc
        h2.get_potential_energy()
        # Context
        calc._pause_calc()
        # Close statement should not last more than 10 sec
        with time_limit(10):
            calc.close()
    return


def test_paused_close_context():
    """Context mode"""
    skip_probe(4)
    skip_slurm()
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(directory=tmpdir, **params) as calc:
            h2.calc = calc
            # Context
            with calc._ensure_mpi():
                h2.get_potential_energy()
            assert calc.process is not None
            pid = calc.process.pid
            assert get_average_cpu(pid) <= 5.0
        # Close statement should not last more than 10 sec
    return


def test_always_true():
    pass
