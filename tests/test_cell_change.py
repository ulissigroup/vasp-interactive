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
    "e": [-3.1472728700846733, -3.0477602406694064, -2.3237350280728015],
    "s": [
        np.array(
            [
                7.46862591e-03,
                7.46846034e-03,
                7.46864299e-03,
                -1.34228406e-05,
                -1.33775906e-05,
                -1.34353895e-05,
            ]
        ),
        np.array(
            [
                6.02763847e-02,
                6.02763500e-02,
                6.02763758e-02,
                -2.90205637e-05,
                -2.90194473e-05,
                -2.90169892e-05,
            ]
        ),
        np.array(
            [
                6.48667580e-02,
                6.48667485e-02,
                6.48667465e-02,
                -3.56225069e-05,
                -3.56320834e-05,
                -3.55993884e-05,
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
            atoms.calc = Vasp(directory="cell-change-vasp", **params)
            e = atoms.get_potential_energy()
            s = atoms.get_stress()
            e_vasp.append(e)
            s_vasp.append(s)
    e_vasp = np.array(e_vasp)
    s_vasp = np.array(s_vasp)

    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(directory="cell-change", **params) as calc:
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
    assert np.linalg.norm(e_vpi - ipi_reference["e"]) < 1.0e-4
    assert np.linalg.norm(s_vpi - ipi_reference["s"]) < 1.0e-5
    # print(e_vasp, s_vasp)
    # print(e_vpi, s_vpi)
