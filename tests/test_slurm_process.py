"""Test process recovery on slurm systems (currenly only NERSC Cori)
"""
import time
import pytest
import numpy as np
from vasp_interactive import VaspInteractive
from vasp_interactive.utils import _locate_slurm_step, _get_slurm_jobid
import psutil
import tempfile
from pathlib import Path
import os
from ase.build import molecule
from ase.optimize import BFGS
from copy import copy, deepcopy
from _common_utils import skip_probe, skip_slurm


h2_root = molecule("H2", vacuum=4, pbc=True)
params = dict(xc="pbe", kpts=(1, 1, 1), ismear=0)
rootdir = Path(__file__).parents[1] / "sandbox"


def test_jobid():
    """Simple test if slurm jobs parsing is ok"""
    skip_slurm(reverse=True)
    h2 = h2_root.copy()
    jobid = _get_slurm_jobid()
    assert jobid is not None
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            h2.calc = calc
            h2.get_potential_energy()
            # Now should have a running VASP srun
            stepid = _locate_slurm_step()
            assert stepid is not None
            assert "." in stepid
            print("jobid", jobid, "stepid", stepid)
            stem, seq = stepid.split(".")
            assert jobid == stem
        # breakpoint()
        # vasp process is terminated
        # wait a few seconds until step finishes
        time.sleep(10)
        stepid = _locate_slurm_step()
        assert stepid is None


def get_average_cpu(pid, interval=0.5):
    proc = psutil.Process(pid)
    vasp_procs = [p for p in proc.children(recursive=True) if "vasp" in p.name()]
    cpu_per = [p.cpu_percent(interval) for p in vasp_procs]
    return np.mean(cpu_per)


def test_signal_send():
    """Simple test if slurm jobs parsing is ok"""
    skip_slurm(reverse=True)
    h2 = h2_root.copy()
    jobid = _get_slurm_jobid()
    assert jobid is not None
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            h2.calc = calc
            h2.get_potential_energy()
            stepid1 = _locate_slurm_step()
            h2.calc._pause_calc()
            time.sleep(10)
            stepid2 = _locate_slurm_step()
            h2.calc._resume_calc()
            stepid3 = _locate_slurm_step()
        assert stepid1 == stepid2 == stepid3
