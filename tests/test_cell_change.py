"""Test energy and stress change of arbitrary structure
"""
import ase
import numpy as np
import tempfile
from ase.build import bulk
from ase.calculators.vasp import Vasp

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible

al = bulk("Al")
# Very bad parameter only for speed
params = {
    "xc": "pbe",
    "encut": 200,
    "istart": 0,
    "lwave": False,
    "ediff": 1.0e-4,
    "kpts": [3, 3, 3],
    "gamma": True,
    "isym": 0,
}
images = []
for r in [1.00, 1.05, 1.20]:
    al_scaled = al.copy()
    al_scaled.set_cell(al.cell * [r, r, r], scale_atoms=True)
    images.append(al_scaled)

ipi_reference = {
    "e": [-3.147273182334173, -3.047760183166516, -2.323741004889209],
    "s": [
        np.array(
            [
                7.47393574e-03,
                7.47235425e-03,
                7.47266351e-03,
                -1.31302173e-05,
                -1.31536271e-05,
                -1.32826985e-05,
            ]
        ),
        np.array(
            [
                6.02688661e-02,
                6.02688401e-02,
                6.02688530e-02,
                -2.90121171e-05,
                -2.90232764e-05,
                -2.90260096e-05,
            ]
        ),
        np.array(
            [
                6.48745402e-02,
                6.48746176e-02,
                6.48745393e-02,
                -3.52549640e-05,
                -3.58915830e-05,
                -3.52846098e-05,
            ]
        ),
    ],
}


def test_cell():
    skip_lattice_if_incompatible()
    e_vasp = []
    s_vasp = []
    e_vpi = []
    s_vpi = []
    # Change lattice scale. Each cell has new vasp calculator
    for atoms in images:
        with tempfile.TemporaryDirectory() as tmpdir:
            atoms.calc = Vasp(directory=tmpdir, **params)
            e = atoms.get_potential_energy()
            s = atoms.get_stress()
            e_vasp.append(e)
            s_vasp.append(s)
    e_vasp = np.array(e_vasp)
    s_vasp = np.array(s_vasp)

    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(directory=tmpdir, **params) as calc:
            for atoms in images:
                atoms.calc = calc
                e = atoms.get_potential_energy()
                s = atoms.get_stress()
                e_vpi.append(e)
                s_vpi.append(s)
    e_vpi = np.array(e_vpi)
    s_vpi = np.array(s_vpi)
    # There might be a tiny difference between single point vasp
    # and vpi, but the difference between vpi and ipi are much smaller
    assert np.linalg.norm(e_vasp - e_vpi) < 0.05
    assert np.linalg.norm(s_vasp - s_vpi) < 5.0e-3
    assert np.linalg.norm(e_vpi - ipi_reference["e"]) < 5.0e-4
    assert np.linalg.norm(s_vpi - ipi_reference["s"]) < 1.0e-4
    # print(e_vasp, s_vasp)
    # print(e_vpi, s_vpi)


def test_always_true():
    pass
