"""Test parser for VASP MLFF
"""
import pytest
from ase.units import GPa
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from vasp_interactive.vasp_interactive import (
    parse_vaspout_energy,
    parse_vaspout_forces,
    _int_version,
)
from vasp_interactive.utils import _preprocess_mlff_outcar
from ase.calculators.calculator import CalculatorSetupError
import tempfile
from pathlib import Path
import os
from ase.io import read
from ase.build import molecule


curdir = Path(__file__).parent
outfile_root = curdir / "mlff_outputs"
atoms = read(outfile_root / "POSCAR")


def test_mlff_flag():
    """Test if VaspInteractive recognizes MLFF flags"""
    vpi = VaspInteractive()
    assert vpi._use_mlff() is False

    # Does not support ML_LMLFF as default params
    with pytest.raises(Exception):
        vpi = VaspInteractive(ml_lmlff=True)

    vpi = VaspInteractive(custom=dict(ml_lmlff=True))
    assert vpi._use_mlff() is True

    vpi = VaspInteractive(custom=dict(ml_lmlff=False))
    assert vpi._use_mlff() is False


def test_ml_istart0():
    """ML_ISTART=0: both DFT calculation and ML inference steps are present
    For the test OUTCAR, the last step is DFT
    """
    refs = {"E": -6.69241300, "F": 0.790714, "S": 0.57098  * -1.e-1 * GPa}  # only z-direction
    vpi = VaspInteractive(custom=dict(ml_lmlff=True))
    # Must use initialize to calculate sort/resort vectors
    vpi.initialize(atoms)
    print(vpi.atoms, vpi.sort, vpi.resort)
    vpi.process = None
    outcar = open(outfile_root / "OUTCAR.istart0", "r").readlines()
    # Last energy is DFT, will work without filtering
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    s = vpi.read_stress(lines=outcar)
    print(s)
    assert e == pytest.approx(refs["E"])
    assert fe == pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])
    # Do pre-filtering
    outcar = _preprocess_mlff_outcar(outcar)
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    assert e == pytest.approx(refs["E"])
    assert fe == pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])



def test_ml_istart1():
    """ML_ISTART=1: both DFT calculation and ML inference steps are present
    For the test OUTCAR, the last step is ML
    """
    refs = {"E": -6.71101275, "F": -4.634836, "S": -10.86173  * -1.e-1 * GPa}  # only z-direction
    vpi = VaspInteractive(custom=dict(ml_lmlff=True))
    # Must use initialize to calculate sort/resort vectors
    vpi.initialize(atoms)
    print(vpi.atoms, vpi.sort, vpi.resort)
    vpi.process = None
    outcar = open(outfile_root / "OUTCAR.istart1", "r").readlines()
    # Last energy is DFT, will work without filtering except for forces
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    s = vpi.read_stress(lines=outcar)
    assert e != pytest.approx(refs["E"])
    assert fe != pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])
    # Do pre-filtering
    outcar = _preprocess_mlff_outcar(outcar)
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    assert e == pytest.approx(refs["E"])
    assert fe == pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])


def test_ml_istart2():
    """ML_ISTART=2: only ML inference
    For the test OUTCAR, the last step is ML
    """
    refs = {"E": -6.75672861, "F": -5.427111, "S": -12.68955  * -1.e-1 * GPa}  # only z-direction
    vpi = VaspInteractive(custom=dict(ml_lmlff=True))
    # Must use initialize to calculate sort/resort vectors
    vpi.initialize(atoms)
    print(vpi.atoms, vpi.sort, vpi.resort)
    vpi.process = None
    outcar = open(outfile_root / "OUTCAR.istart2", "r").readlines()
    # Last energy is DFT, will work without filtering except for forces
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    s = vpi.read_stress(lines=outcar)
    assert e != pytest.approx(refs["E"])
    assert fe != pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])
    # Do pre-filtering
    outcar = _preprocess_mlff_outcar(outcar)
    e, fe = vpi.read_energy(lines=outcar)
    f = vpi.read_forces(lines=outcar)
    assert e == pytest.approx(refs["E"])
    assert fe == pytest.approx(refs["E"])
    assert f[0, 2] == pytest.approx(refs["F"])
    assert s[2] == pytest.approx(refs["S"])
