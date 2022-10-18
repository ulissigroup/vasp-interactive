"""Compare the energy and stress outputs from VASP IBRION+ISIF relaxation
with VaspInteractive's stress mode
"""
import ase
import numpy as np
import tempfile
import pickle
import pytest
from pathlib import Path
from ase.io import read, write
from ase.build import molecule
from ase.calculators.vasp import Vasp
from ase.calculators.socketio import SocketServer, SocketIOCalculator

from vasp_interactive import VaspInteractive
from _common_utils import skip_lattice_if_incompatible

water = molecule("H2O", vacuum=5, pbc=True)
params = {"xc": "pbe", "encut": 150, "ediff": 1.e-3, "istart": 0, "lwave": False}

def test_socket_dryrun():
    """Initialize the VaspInteractive calculator without running socket
    """
    # Case 1: client without initialization
    vpi = VaspInteractive(**params)
    water.calc = vpi
    print(vpi.command)
    assert "-m vasp_interactive.socketio" in vpi.command
    assert "-ht localhost" in vpi.command
    assert vpi.socket_client is None
    # Does the vpi initial parameter exists?
    param_file = Path(vpi._indir(".vpi_params.pkl"))
    assert param_file.is_file()
    input_params = pickle.load(open(param_file, "rb"))
    assert all([p not in input_params for p in ("use_socket", "port", "unixsocket", "host")])

    # Case 2: add socket without server, raises error
    with pytest.raises(ConnectionRefusedError):
        vpi = VaspInteractive(use_socket=True, **params)
    
    # Case 3: create server, do not run
    server = SocketServer(port=31415)
    # client should connect without issue
    vpi = VaspInteractive(use_socket=True, port=31415, **params)
    # terminate the client and server
    assert vpi.socket_client is not None
    assert vpi.socket_client.state == "READY"
    vpi.finalize()
    server.close()
    assert vpi.socket_client.closed
    assert vpi.final

def test_socket_port():
    vpi = VaspInteractive(**params)
    server = SocketIOCalculator(calc=vpi, port=31415)
    # Using the wrapper, vpi is not directly used as a calculator
    assert server.server is None
    assert vpi.socket_client is None
    atoms = water.copy()
    atoms.calc = server
    with server:
        for i in range(2):
            atoms.rattle()
            # server.calculate(atoms)
            e = atoms.get_potential_energy()
            assert server.server is not None
            assert server.server.proc is not None
            assert vpi.socket_client is None
            assert vpi.process is None
            print(e)
            # print(res)

def test_socket_unixsocket():
    vpi = VaspInteractive(**params)
    server = SocketIOCalculator(calc=vpi, unixsocket="vpi")
    # Using the wrapper, vpi is not directly used as a calculator
    assert server.server is None
    assert vpi.socket_client is None
    atoms = water.copy()
    atoms.calc = server
    with server:
        for i in range(2):
            atoms.rattle()
            # server.calculate(atoms)
            e = atoms.get_potential_energy()
            assert server.server is not None
            assert server.server.proc is not None
            assert vpi.socket_client is None
            assert vpi.process is None
            print(e)