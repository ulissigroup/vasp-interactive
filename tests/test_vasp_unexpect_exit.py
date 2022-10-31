"""Test if VaspInteractive can handle unexpected VASP close
"""
import pytest
import numpy as np
from vasp_interactive import VaspInteractive
import psutil
import time
import tempfile
from subprocess import run, Popen
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



def test_abrupt_kill():
    """Randomly kill vasp process during a run
    """
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, txt="-", **params)
        with calc:
            sleep_t = 3 + np.random.random() * 2
            # Test if abruptly killing vasp causes RuntimeError
            p = Popen(f"sleep {sleep_t} && killall vasp_std", shell=True)
            h2.calc = calc
            try:
                h2.get_potential_energy()
            except RuntimeError:
                assert p.poll() is not None
    return


# def test_pause_context():
#     """Context mode"""
#     skip_probe(4, skip_oversubscribe=True)
#     # skip_slurm()
#     h2 = h2_root.copy()
#     threshold_high = 75.0
#     threshold_low = 25.0
#     with tempfile.TemporaryDirectory() as tmpdir:
#         with VaspInteractive(directory=tmpdir, **params) as calc:
#             h2.calc = calc
#             h2.get_potential_energy()
#             pid = h2.calc.process.pid
#             # Context
#             with h2.calc.pause():
#                 cpu_stop = get_average_cpu()
#                 print(cpu_stop)
#                 assert cpu_stop < threshold_low

#             cpu_nonstop = get_average_cpu()
#             print(cpu_nonstop)
#             assert cpu_nonstop > threshold_high
#     return


# def test_always_true():
#     """Placeholder"""
#     pass
