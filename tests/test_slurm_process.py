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
from _common_utils import skip_probe, skip_slurm, get_average_cpu


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
            pid1 = calc.process.pid
            stepid1 = _locate_slurm_step()
            h2.calc._pause_calc()
            time.sleep(2)
            cpu = get_average_cpu()
            print(f"After pause, cpu is {cpu}")
            assert cpu < 25
            pid2 = calc.process.pid
            stepid2 = _locate_slurm_step()
            h2.calc._resume_calc()
            cpu = get_average_cpu()
            assert cpu > 75
            print(f"After cont, cpu is {cpu}")
            pid3 = calc.process.pid
            stepid3 = _locate_slurm_step()
        assert stepid1 == stepid2 == stepid3
        assert pid1 == pid2 == pid3
        with calc:
            h2.rattle(0.1)
            h2.get_potential_energy()
            pid4 = calc.process.pid
            stepid4 = _locate_slurm_step()
            assert pid1 != pid4
            assert stepid1 != stepid4


def test_signal_time():
    """Timing test for slurm pause / resume
    By design the second pause / resume step should take much less time than first one
    """
    skip_slurm(reverse=True)
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(directory=tmpdir, **params)
        with calc:
            for i in range(3):
                h2.calc = calc
                h2.rattle(0.1)
                h2.get_potential_energy()
                if i == 0:
                    assert getattr(h2.calc, "mpi_match", None) is None
                else:
                    assert getattr(h2.calc, "mpi_match", None) is not None
                t_ = time.time()
                h2.calc._pause_calc()
                assert getattr(h2.calc, "mpi_match", None) is not None
                t_el1 = time.time() - t_
                print(f"Pause {i}, time {t_el1} s")
                time.sleep(2)
                t_ = time.time()
                h2.calc._resume_calc()
                t_el2 = time.time() - t_
                print(f"Cont {i}, time {t_el2} s")
                # Loose condition to allow unexpected slurm turnaround degredation
                # 30 sec timeout should normally be met for querying `squeue`
                if i == 0:
                    assert t_el1 < 30
                    assert t_el2 < 10
                else:
                    assert t_el1 < 10
                    assert t_el2 < 10
