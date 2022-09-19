"""Use VaspInteractive to calculate relaxation of H2 molecule
"""

import os
import sys
import subprocess
import tempfile
from ase.build import molecule
from ase.io import read
from vasp_interactive import VaspInteractive


def vasp_env_test():
    print(
        "Checking if ASE_VASP_COMMAND, VASP_COMMAND or VASP_SCRIPT environment variables are set"
    )
    if any(
        [
            os.environ.get(p, None) is not None
            for p in ("ASE_VASP_COMMAND", "VASP_COMMAND", "VASP_SCRIPT")
        ]
    ):
        print("env set OK!")
    else:
        var = input("Please enter the VASP command:\n")
        if len(var) == 0:
            print("Failed to set VASP env. Abort.")
            return False
        else:
            os.environ["ASE_VASP_COMMAND"] = str(var)
            print(f"set ASE_VASP_COMMAND to {var}")
    print("Checking VASP Pseudopotential settings")

    if os.environ.get("VASP_PP_PATH", None) is not None:
        print("VASP_PP_PATH set")
    else:
        var = input("Please enter the VASP_PP_PATH environmen variable:\n")
        if len(var) == 0:
            print("Failed to set VASP_PP_PATH env. Abort.")
            return False
        else:
            os.environ["VASP_PP_PATH"] = str(var)
            print(f"set VASP_PP_PATH to {var}")
    return True


def demo_test():
    print("Running test example")
    atoms = molecule("H2", pbc=True, vacuum=4)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            xc="pbe",
            directory="test1",
        )
        # Low level calculator interfacing
        with calc._txt_outstream() as out:
            calc._run(atoms, out=out)
        pid = calc.process.pid

        # Check vasprun.xml

        # Check OUTCAR
        outcar_lines = open(calc._indir("OUTCAR"), "r").readlines()
        outcar_ok = False
        cond = 0
        matches = [
            "FORCE on cell",
            "TOTAL-FORCE",
            "energy  without entropy",
            "VOLUME and BASIS-vectors",
            "E-fermi",
        ]
        for line in outcar_lines:
            if any([m in line for m in matches]):
                cond += 1
        if cond >= 5:
            outcar_ok = True

        # Check vaspout
        vaspout_lines = calc._txt_to_handler().readlines()
        vaspout_ok = False
        for line in vaspout_lines:
            if "FORCES" in line:
                vaspout_ok = True
                break
        # Low level kill
        subprocess.run(["kill", "-9", str(pid)])

    print(f"VASP OUTCAR: {outcar_ok}")
    print(f"VASP raw output: {vaspout_ok}")


if __name__ == "__main__":
    if not vasp_env_test():
        sys.exit(1)
    demo_test()
