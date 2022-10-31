import sys

from ase.build import molecule, bulk
from ase.io import write
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
from subprocess import Popen
import os

# os.environ["VASP_COMMAND"] = "mpirun -np 8 --map-by hwthread /home/jovyan/data/vasp_binaries/patch/vasp63/bin/vasp_std"
os.environ[
    "VASP_COMMAND"
] = "mpirun -np 8 --map-by hwthread /home/jovyan/data/vasp_binaries/pristine/vasp63/bin/vasp_std"

al = bulk("Al")
al_scaled = al.copy()
# Very bad parameter only for speed
params = {
    "xc": "pbe",
    "encut": 200,
    "ediff": 1.0e-4,
    "kpts": [3, 3, 3],
    "gamma": True,
    "lwave": False,
    "istart": 0,
}
#  "custom": {"ihost": "localhost", "port": "23333"}}


unixsocket = "localhost"

# atoms = molecule('H2O', vacuum=3.0)
# atoms.rattle(stdev=0.1)
# write('initial.traj', atoms)


# with SocketIOCalculator(log=sys.stdout, port=23333) as calc:
# with SocketIOCalculator(log=sys.stdout, unixsocket=unixsocket) as calc:
# Server is now running and waiting for connections.
# If you want to launch the client process here directly,
# instead of manually in the terminal, uncomment these lines:
#
# from subprocess import Popen
# proc = Popen(os.environ["VASP_COMMAND"], cwd="inet-test", shell=True)
vp = VaspInteractive(
    directory="inet-test", use_socket=True, unixsocket="localhost", **params
)
# al_scaled.calc = vp
vp.run(al_scaled)
# print("here")
# blocking
# vp.run(al_scaled)
# print("here2")

#     for r in [1.00, 1.05, 1.20]:
#         al_scaled.set_cell(al.cell * [r, r, r], scale_atoms=True)
#         e = al_scaled.get_potential_energy()
#         s = al_scaled.get_stress()
#         es.append(e)
#         ss.append(s)
# print(es)
# print(ss)
# print(r)
# print(e)
# print(s)
