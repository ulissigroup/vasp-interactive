import pytest
import psutil
import tempfile
import numpy as np
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from ase.atoms import Atoms
import re


def get_cpu_cores():
    cores = 1
    with VaspInteractive() as test_calc:
        args = test_calc.make_command().split()
        for i, param in enumerate(args):
            if param in ["-n", "-np", "--n", "--np", "-c"]:
                try:
                    cores = int(args[i + 1])
                    break
                except Exception as e:
                    cores = 1
                    break
    return cores


def get_oversubscribe():
    over = False
    with VaspInteractive() as test_calc:
        args = test_calc.make_command().split()
        for i, param in enumerate(args):
            if "oversubscribe" in param:
                over = True
                break
    return over


def if_vasp_accepts_lattice():
    """Run a simple test to see if VASP accepts lattice input"""
    h = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(
            directory=tmpdir, xc="pbe", encut=150, ediff=1.0e-3
        ) as calc:
            h.calc = calc
            for i in range(2):
                h.rattle(0.02)
                h.get_potential_energy()
            f_vaspout = calc._txt_to_handler()
            lines = f_vaspout.readlines()
        accept = False
        for line in lines:
            if "LATTICE: reading from stdin" in line:
                accept = True
                break
    return accept


def if_vasp_ipi():
    """Test if vasp has iPi patch by parsing the stderr"""
    h = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(
            directory="./test-ipi",
            xc="pbe",
            ibrion=23,
            custom={"ihost": "localhost", "port": "23333", "inet": "1"},
        )
        cmd = calc.make_command(calc.command) + " 2> vasp.err"
        calc.command = cmd
        h.calc = calc
        try:
            h.get_potential_energy()
        except Exception as e:
            pass
        outcar_lines = open(calc._indir("vasp.err"), "r").readlines()
        accept = False
        for line in outcar_lines:
            if "Error opening INET socket" in line:
                accept = True
                break
    return accept


def skip_probe(min_cores=8, skip_oversubscribe=False):
    """Test if single step needs to be skipped"""
    cores = get_cpu_cores()
    do_test = cores >= min_cores
    if do_test is False:
        pytest.skip(
            f"Skipping test with ncores < {min_cores}", allow_module_level=False
        )
    elif skip_oversubscribe and get_oversubscribe():
        pytest.skip(f"Skipping due to oversubscription", allow_module_level=False)


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


def skip_lattice_if_incompatible(reverse=False):
    """Skip tests that requires cell changes
    the tests automatically detect if VASP is compatible with lattice input
    """
    accept = if_vasp_accepts_lattice()
    # skip if accept and reverse are same logic value
    if not (accept ^ reverse):
        pytest.skip(
            f"Skipping because vasp does not supports lattice input",
            allow_module_level=False,
        )


def skip_if_no_ipi(reverse=False):
    """Skip tests that requires ipi patch"""
    accept = if_vasp_ipi()
    # skip if accept and reverse are same logic value
    if not (accept ^ reverse):
        pytest.skip(
            f"Skipping because vasp does not support ipi protocol",
            allow_module_level=False,
        )
