import pytest
import psutil
import numpy as np
from vasp_interactive import VaspInteractive


def skip_probe(min_cores=8):
    """Test if single step needs to be skipped"""
    with VaspInteractive() as test_calc:
        args = test_calc.make_command().split()
        do_test = True
        for i, param in enumerate(args):
            if param in ["-n", "-np", "--n", "--np", "-c"]:
                try:
                    cores = int(args[i + 1])
                    do_test = cores >= min_cores
                    break
                except Exception as e:
                    do_test = False
                    break
        if do_test is False:
            pytest.skip(
                f"Skipping test with ncores < {min_cores}", allow_module_level=False
            )


def skip_slurm(reverse=False):
    """Test if single step needs to be skipped"""
    with VaspInteractive() as test_calc:
        args = test_calc.make_command().split()
        do_test = True
        for i, param in enumerate(args):
            if "srun" in param:
                do_test = False
                break
    if reverse:
        do_test = not do_test
    if do_test is False:
        pytest.skip(f"Skipping test started by srun", allow_module_level=False)


def get_average_cpu(interval=0.5):
    """Get average cpu usage for vasp processes"""
    vasp_procs = [p for p in psutil.process_iter() if "vasp" in p.name()]
    cpu_per = [p.cpu_percent(interval) for p in vasp_procs]
    return np.mean(cpu_per)
