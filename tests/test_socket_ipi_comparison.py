"""Compare the results of E, F and S between VaspInteractive socket IO and 
iPI-patch
"""
import ase
import numpy as np
import tempfile
import sys
import pickle
import pytest
from multiprocessing import Process
from pathlib import Path
from ase.io import read, write
from ase.build import molecule
from ase.calculators.vasp import Vasp
from ase.calculators.socketio import SocketServer, SocketIOCalculator

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible, skip_if_no_ipi

curdir = Path(__file__).parent.resolve()
water = molecule("H2O", vacuum=5, pbc=True)
gaas = read(curdir / "structures" / "POSCAR.GaAs", format="vasp")

params_water = {"xc": "pbe", "encut": 150, "ediff": 1.0e-3, "istart": 0, "lwave": False}
params_gaas = {
    "xc": "pbe",
    "encut": 300,
    "ediff": 1.0e-4,
    "istart": 0,
    "kpts": [6, 6, 3],
    "sigma": 0.1,
}


def test_molecule():
    """Run socket client in a separate Process at background"""
    skip_if_no_ipi()
    e_vpi, e_ipi = [], []
    f_vpi, f_ipi = [], []
    n_samples = 3
    images = []
    for i in range(n_samples):
        atoms = water.copy()
        atoms.rattle(0.02)
        images.append(atoms.copy())
    # 1. VaspInteractive
    with tempfile.TemporaryDirectory() as tempdir:
        vpi = VaspInteractive(directory=tempdir, **params_water)
        server = SocketIOCalculator(calc=vpi, log=sys.stdout, unixsocket="vpi")
        # Using the wrapper, vpi is not directly used as a calculator
        with server:
            for i in range(n_samples):
                atoms = images[i]
                atoms.calc = server
                # server.calculate(atoms)
                e = atoms.get_potential_energy()
                f = atoms.get_forces()
                e_vpi.append(e)
                f_vpi.append(f)
    # 2. iPI-VASP
    with tempfile.TemporaryDirectory() as tempdir:
        vpi = Vasp(
            directory=tempdir,
            custom={"port": "2333", "inet": "1"},
            ibrion=23,
            **params_water
        )
        vpi.command = vpi.make_command()
        server = SocketIOCalculator(calc=vpi, log=sys.stdout, port=2333)
        # Using the wrapper, vpi is not directly used as a calculator
        with server:
            for i in range(n_samples):
                atoms = images[i]
                atoms.calc = server
                e = atoms.get_potential_energy()
                f = atoms.get_forces()
                e_ipi.append(e)
                f_ipi.append(f)
    # print(e_vpi, e_ipi)
    print(f_vpi)
    print(f_ipi)
    assert np.linalg.norm(np.array(e_vpi) - np.array(e_ipi)) < 1.0e-3
    assert np.linalg.norm(np.array(f_vpi) - np.array(f_ipi), axis=-1).max() < 1.0e-3


def test_bulk():
    """Run socket client in a separate Process at background"""
    skip_if_no_ipi()
    skip_lattice_if_incompatible()
    e_vpi, e_ipi = [], []
    f_vpi, f_ipi = [], []
    s_vpi, s_ipi = [], []
    n_samples = 3
    images = []
    np.random.seed(42)
    for i in range(n_samples):
        atoms = gaas.copy()
        atoms.rattle(0.02)
        strain = 0.1 * (np.random.random() - 0.5)
        cell = atoms.cell * [1 + strain, 1 + strain, 1 + strain]
        atoms.set_cell(cell, scale_atoms=True)
        images.append(atoms.copy())
    # 1. VaspInteractive
    with tempfile.TemporaryDirectory() as tempdir:
        vpi = VaspInteractive(directory=tempdir, **params_water)
        server = SocketIOCalculator(calc=vpi, log=sys.stdout, unixsocket="vpi")
        # Using the wrapper, vpi is not directly used as a calculator
        with server:
            for i in range(n_samples):
                atoms = images[i]
                atoms.calc = server
                # server.calculate(atoms)
                e = atoms.get_potential_energy()
                f = atoms.get_forces()
                s = atoms.get_stress()
                e_vpi.append(e)
                f_vpi.append(f)
                s_vpi.append(s)
    # 2. iPI-VASP
    with tempfile.TemporaryDirectory() as tempdir:
        vpi = Vasp(
            directory=tempdir,
            custom={"port": "2333", "inet": "1"},
            ibrion=23,
            **params_water
        )
        vpi.command = vpi.make_command()
        server = SocketIOCalculator(calc=vpi, log=sys.stdout, port=2333)
        # Using the wrapper, vpi is not directly used as a calculator
        with server:
            for i in range(n_samples):
                atoms = images[i]
                atoms.calc = server
                e = atoms.get_potential_energy()
                f = atoms.get_forces()
                s = atoms.get_stress()
                e_ipi.append(e)
                f_ipi.append(f)
                s_ipi.append(s)
    # print(e_vpi, e_ipi)
    # print(f_vpi)
    # print(f_ipi)
    # print(s_vpi)
    print(np.array(f_vpi) - np.array(f_ipi))
    print(np.array(s_vpi) - np.array(s_ipi))
    assert np.linalg.norm(np.array(e_vpi) - np.array(e_ipi)) < 1.0e-3
    assert np.linalg.norm(np.array(f_vpi) - np.array(f_ipi), axis=-1).max() < 1.0e-2
    assert np.linalg.norm(np.array(s_vpi) - np.array(s_ipi)) < 1.0e-3
