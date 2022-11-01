"""Test if VaspInteractive can handle unexpected VASP close
"""
import pytest
import numpy as np
from vasp_interactive import VaspInteractive
import psutil
import time
import tempfile
import signal
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


# def test_abrupt_kill():
#     """Randomly kill vasp process during a run"""
#     h2 = h2_root.copy()
#     with tempfile.TemporaryDirectory() as tmpdir:
#         calc = VaspInteractive(directory=tmpdir, txt="-", **params)
#         with calc:
#             sleep_t = np.random.uniform(2, 3)
#             # Test if abruptly killing vasp causes RuntimeError
#             p = Popen(f"sleep {sleep_t} && killall vasp_std", shell=True)
#             h2.calc = calc
#             try:
#                 h2.get_potential_energy()
#             except RuntimeError:
#                 assert p.poll() is not None
#             # Do not exit the function while subprocess is still running
#             while p.poll() is not None:
#                 time.sleep(0.5)
#     return


def test_pause_kill():
    """Randomly kill vasp process during a run"""
    h2 = h2_root.copy()
    tmpdir = tempfile.mkdtemp()
    # with tempfile.TemporaryDirectory() as tmpdir:
    calc = VaspInteractive(directory=tmpdir, txt="-", **params)
    with calc:
        print(calc.pause_mpi)
        h2.calc = calc
        h2.get_potential_energy()
        calc._pause_calc()
        # Manually kill the VASP process
        calc.process.kill()
        # calc._resume_calc should still work
        calc._resume_calc()
        h2.rattle(0.05)
        print(calc.process)
        print(calc.pid)
        # At this stage when calc.process checks it will return None zero exit code
        with pytest.raises(RuntimeError):
            h2.get_potential_energy()
        rt = calc.process.poll()
        assert rt and (rt != 0)
    return


def test_pause_low_nsw():
    """If VASP exits due to a low nsw setting, _resume_calc and _pause_calc
    should not causing errors, but rather start new VASP instance.
    """
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, nsw=1, txt="-", **params)
        with calc:
            h2.calc = calc
            calc._pause_calc()
            calc._resume_calc()
            h2.get_potential_energy()
            pid1 = calc.pid
            # At this step calc already quits.
            calc._pause_calc()
            assert calc.process is None
            assert calc.pid is None
            calc._resume_calc()
            h2.rattle(0.05)
            h2.get_potential_energy()
            pid2 = calc.pid
            assert pid1 != pid2
    return


def test_abrupt_stopcar():
    """Randomly write STOPCAR to stop a VASP calculation.
    pause / resume should not be affected.
    """
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, txt="-", **params)
        with calc:
            sleep_t = np.random.uniform(1, 2)
            # First ionic step only has 1 scf
            h2.calc = calc
            p = Popen(
                f"sleep {sleep_t} && echo 'LABORT = .TRUE.' > STOPCAR",
                shell=True,
                cwd=tmpdir,
            )
            e1 = h2.get_potential_energy()
            print(e1)
            calc._pause_calc()
            calc._resume_calc()
            #
            h2.rattle(0.05)
            e2 = h2.get_potential_energy()
            print(e2)
            print(calc.process)
            # assert calc.process.poll() == 0
            calc._pause_calc()
            calc._resume_calc()
            h2.rattle(0.05)
            e3 = h2.get_potential_energy()
            print(e3)
            # assert e1 == e2 != e3
    return


# def test_always_true():
#     """Placeholder"""
#     pass
