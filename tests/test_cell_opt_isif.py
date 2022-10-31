"""Compare the energy and stress outputs from VASP IBRION+ISIF relaxation
with VaspInteractive's stress mode
"""
import ase
import numpy as np
import tempfile
from ase.io import read, write
from ase.build import bulk
from ase.calculators.vasp import Vasp

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible

a_ = 2.95
al = bulk("Al", a=a_, b=a_, c=a_)
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

with tempfile.TemporaryDirectory() as tmpdir:
    atoms = al.copy()
    atoms.calc = Vasp(directory=tmpdir, isif=3, ibrion=1, nsw=100, **params)
    # Do a direct cell relaxation
    atoms.get_potential_energy()
    outcar = atoms.calc._indir("OUTCAR")
    images = read(outcar, ":")

vasp_es = []
vasp_ss = []
for im in images:
    vasp_es.append(im.get_potential_energy())
    vasp_ss.append(np.max(np.abs(im.get_stress())))
    # print(im.get_potential_energy())
    # print(im.get_stress())
vasp_es = np.array(vasp_es)
vasp_ss = np.array(vasp_ss)
print(vasp_es)
print(vasp_ss)


def test_reference():
    skip_lattice_if_incompatible()
    e_vpi = []
    s_vpi = []

    with tempfile.TemporaryDirectory() as tmpdir:
        with VaspInteractive(directory=tmpdir, **params) as calc:
            for im in images:
                atoms = im.copy()
                atoms.calc = calc
                e = atoms.get_potential_energy()
                s = atoms.get_stress()
                e_vpi.append(e)
                s_vpi.append(np.max(np.abs(s)))
    e_vpi = np.array(e_vpi)
    s_vpi = np.array(s_vpi)
    print(e_vpi)
    print(s_vpi)
    print((e_vpi - vasp_es))
    print((s_vpi - vasp_ss))
    assert np.linalg.norm(e_vpi - vasp_es) < 1.0e-5
    assert np.linalg.norm(s_vpi - vasp_ss) < 1.0e-4


def test_always_true():
    pass
