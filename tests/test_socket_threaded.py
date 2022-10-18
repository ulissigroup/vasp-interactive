"""Compare the energy and stress outputs from VASP IBRION+ISIF relaxation
with VaspInteractive's stress mode
"""
import ase
import numpy as np
import tempfile
import pickle
import pytest
from multiprocessing import Process
from pathlib import Path
from ase.io import read, write
from ase.build import molecule
from ase.calculators.vasp import Vasp
from ase.calculators.socketio import SocketServer, SocketIOCalculator

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible

water = molecule("H2O", vacuum=5, pbc=True)
params = {"xc": "pbe", "encut": 150, "ediff": 1.e-3, "istart": 0, "lwave": False}


def test_socket_unixsocket():
    """Run socket client in a separate Process at background
    """
    server = SocketIOCalculator(calc=None, unixsocket="vpi")
    vpi = VaspInteractive(use_socket=True, unixsocket="vpi", **params)
    assert vpi.socket_client is not None
    assert vpi.socket_client.parent_calc is vpi
    # Using the wrapper, vpi is not directly used as a calculator
    atoms = water.copy()
    atoms.calc = server
    thread = Process(target=vpi.run, args=(atoms,))
    thread.start()
    with server:
        for i in range(2):
            print(i)
            atoms.rattle()
            # server.calculate(atoms)
            e = atoms.get_potential_energy()
            assert vpi.process is None
            # print(e)
    thread.join()
    assert vpi.process is None
    vpi.finalize()

def test_socket_port():
    """Run socket client in a separate Process at background
    """
    server = SocketIOCalculator(calc=None, unixsocket=31415)
    vpi = VaspInteractive(use_socket=True, unixsocket=31415, **params)
    # Using the wrapper, vpi is not directly used as a calculator
    atoms = water.copy()
    atoms.calc = server
    thread = Process(target=vpi.run, args=(atoms,))
    thread.start()
    with server:
        for i in range(2):
            print(i)
            atoms.rattle()
            e = atoms.get_potential_energy()
            assert vpi.process is None
    thread.join()
    assert vpi.process is None
    vpi.finalize()
