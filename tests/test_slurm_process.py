"""Test process recovery on slurm systems (currenly only NERSC Cori)
"""
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
        calc = VaspInteractive(directory=rootdir / "slurm-test", **params)
        with calc:
            h2.calc = calc
            h2.get_potential_energy()
            # Now should have a running VASP srun
            stepid = _locate_slurm_step()
            assert stepid is not None
            print("jobid", jobid, "stepid", stepid)
            assert jobid in stepid
