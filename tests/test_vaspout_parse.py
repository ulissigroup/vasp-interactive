"""Testing if vasp.out parsing works for VASP 5.x
broken_output directory has 2 examples of truncated output blocks from VASP 5.x
step-0
"""
import pytest
from vasp_interactive import VaspInteractive
from vasp_interactive.vasp_interactive import (
    parse_vaspout_energy,
    parse_vaspout_forces,
    _int_version,
)
from ase.calculators.calculator import CalculatorSetupError
import tempfile
from pathlib import Path
import os
from ase.io import read
from ase.build import molecule

curdir = Path(__file__).parent
outfile_root = curdir / "broken_output"


def test_vaspout_energy():
    for case in ["step-0", "step-1"]:
        lines = open(outfile_root / case / "vasp.out", "r").readlines()
        fe, e0 = parse_vaspout_energy(lines)
        assert fe != 0
        assert e0 != 0
        assert fe == pytest.approx(-6.2, 0.1)
        assert e0 == pytest.approx(-6.2, 0.5)
        assert fe == pytest.approx(e0, 1.0e-4)


def test_vaspout_forces():
    for case in ["step-0", "step-1"]:
        lines = open(outfile_root / case / "vasp.out", "r").readlines()
        f = parse_vaspout_forces(lines)
        assert f.shape[0] == 2
        assert f.shape[1] == 3
        # fx[0] and fx[1] should compensate
        assert f[0, 0] + f[1, 0] == pytest.approx(0, 1.0e-6)


def test_outcar_failure():
    """calc.read_energy and calc.read_forces will not work unless vasp5=True provided
    In the case of multiple steps, direct reading from OUTCAR will even mess up the values
    """
    reference = {
        "energy": {"step-0": -0.62323928e01, "step-1": -0.62232316e01},
        "forces": {"step-0": 3.6669864, "step-1": 3.6854838},
    }
    for case in ["step-0", "step-1"]:
        dir = outfile_root / case
        calc = VaspInteractive(directory=dir)
        # No vasp5, fail
        fe, e0 = calc.read_energy(all=False, vasp5=False)
        assert fe != reference["energy"][case]
        f = calc.read_forces(all=False, vasp5=False)
        if f is not None:
            assert abs(f[0, 0]) != reference["forces"][case]
        # vasp5=True, ok
        fe, e0 = calc.read_energy(all=False, vasp5=True)
        assert fe == pytest.approx(reference["energy"][case], 1.0e-6)
        f = calc.read_forces(all=False, vasp5=True)
        assert f is not None
        assert abs(f[0, 0]) == pytest.approx(reference["forces"][case], 1.0e-6)


def test_calc_option():
    """Higher level test for parse_vaspout option"""
    reference = {
        "energy": {"step-0": -0.62323928e01, "step-1": -0.62232316e01},
        "forces": {"step-0": 3.6669864, "step-1": 3.6854838},
    }
    # No parsing of vasp.out, fail
    for case in ["step-0", "step-1"]:
        dir = outfile_root / case
        calc = VaspInteractive(directory=dir, parse_vaspout=False)
        calc.read_results()
        fe = calc.results["free_energy"]
        f = calc.results["forces"]
        assert fe != reference["energy"][case]
        if f is not None:
            assert abs(f[0, 0]) != reference["forces"][case]

    for case in ["step-0", "step-1"]:
        dir = outfile_root / case
        calc = VaspInteractive(directory=dir, parse_vaspout=True)
        calc.read_results()
        fe = calc.results["free_energy"]
        f = calc.results["forces"]
        assert fe == pytest.approx(reference["energy"][case], 1.0e-6)
        assert f is not None
        assert abs(f[0, 0]) == pytest.approx(reference["forces"][case], 1.0e-6)


def test_txt_option():
    """Parsing with vasp 5 should only work when the txt output of VASP calculator is not suppressed"""
    atoms = molecule("H2", vacuum=5, pbc=True)
    # Very rough settings
    params = dict(
        xc="pbe", kpts=(1, 1, 1), nsw=0, ibrion=-1, ismear=0, ediff=1e-2, encut=120
    )
    # Change vasp.out name
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms0 = atoms.copy()
        calc = VaspInteractive(directory=tmpdir, txt="vasp.out.other", **params)
        atoms0.calc = calc
        fe = atoms0.get_potential_energy()
        version = _int_version(calc.version)
        if version == 6:
            pytest.skip(f"Skipping vasp 6 tests", allow_module_level=False)

        assert fe == pytest.approx(-6.2, 0.1)

    # stdout
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms0 = atoms.copy()
        calc = VaspInteractive(directory=tmpdir, txt="-", **params)
        atoms0.calc = calc
        with pytest.raises(RuntimeError):
            fe = atoms0.get_potential_energy()

    # supressed
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms0 = atoms.copy()
        calc = VaspInteractive(directory=tmpdir, txt=None, **params)
        atoms0.calc = calc
        with pytest.raises(RuntimeError):
            fe = atoms0.get_potential_energy()
