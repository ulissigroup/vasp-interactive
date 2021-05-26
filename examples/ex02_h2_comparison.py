"""Compare results for relaxation of H2 molecule
   using VaspInteractive vs Vasp internal vs Vasp SinglePoint + BFGS
   using the same force stop criteria the energy results should be almost identical
"""
import numpy as np
import os
import tempfile

from ase.atoms import Atoms
from ase.optimize import BFGS

from vasp_interactive import VaspInteractive
from ase.calculators.vasp import Vasp

d = 0.9575
h2_root = Atoms("H2", 
           positions=[(d, 0, 0), (0, 0, 0)], 
           cell=[8, 8, 8])

def run_with_vasp_interactive():
    h2 = h2_root.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(ismear=0,
                               xc="pbe",
                               kpts=(1, 1, 1),  # not important, just keeps it faster
                               directory=tmpdir
                              )
      
        # Best practice of VaspInteractive is to use it as ContextManager
        with calc:
            h2.set_calculator(calc)
            dyn = BFGS(h2)
            # Now ASE-BFGS controls the relaxation, not VASP
            dyn.run(fmax=0.05)
            print("Final energy using VaspInteractive: ", h2.get_potential_energy())

def run_with_vasp():
    h2 = h2_root.copy()
    h2.pbc = True
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(ismear=0,
                    xc="pbe",
                    kpts=(1, 1, 1),  # not important, just keeps it faster
                    directory=tmpdir,
                    ibrion=2,
                    nsw=10,
                    ediffg=-0.05,
                   )
        # Best practice of VaspInteractive is to use it as ContextManager
        h2.set_calculator(calc)
        print("Final energy from Vasp internal routines: ", h2.get_potential_energy())
        
def run_with_vasp_bfgs():
    h2 = h2_root.copy()
    h2.pbc = True
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = Vasp(ismear=0,
                    xc="pbe",
                    kpts=(1, 1, 1),  # not important, just keeps it faster
                    directory=tmpdir,
                    ibrion=-1,
                    nsw=0,)
        # Best practice of VaspInteractive is to use it as ContextManager
        h2.set_calculator(calc)
        dyn = BFGS(h2)
        # Now ASE-BFGS controls the relaxation, not VASP
        dyn.run(fmax=0.05)
        print("Final energy from Vasp + BFGS: ", h2.get_potential_energy())

if __name__ == "__main__":
    run_with_vasp()
    run_with_vasp_bfgs()
    run_with_vasp_interactive()