"""Use VaspInteractive to calculate relaxation of H2 molecule
"""

import os
import sys
import subprocess
import tempfile
from ase.build import molecule
from ase.io import read
from vasp_interactive import VaspInteractive
from vasp_interactive.utils import time_limit
import warnings
import traceback


warnings.filterwarnings("ignore")


def cprint(content, color=None, **kwargs):
    """Color print wrapper"""
    ansi_color = dict(
        HEADER="\033[95m",
        OKBLUE="\033[94m",
        OKGREEN="\033[92m",
        WARNING="\033[93m",
        FAIL="\033[91m",
        ENDC="\033[0m",
        BOLD="\033[1m",
        UNDERLINE="\033[4m",
    )
    if color is None:
        output = content
    elif color in ansi_color.keys() and color != "ENDC":
        output = ansi_color[color] + content + ansi_color["ENDC"]
    else:
        raise ValueError(
            f"Unknown ANSI color name. Allowed values are {list(ansi_color.keys())}"
        )
    print(output, **kwargs)
    return


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
        cprint("VASP command env is set!", color="OKGREEN")
    else:
        var = input("Please enter the VASP command:\n")
        if len(var) == 0:
            cprint("Failed to set VASP env. Abort.", color="FAIL")
            return False
        else:
            os.environ["ASE_VASP_COMMAND"] = str(var)
            print(f"set ASE_VASP_COMMAND to {var}")
    print("Checking VASP Pseudopotential settings")

    if os.environ.get("VASP_PP_PATH", None) is not None:
        cprint("VASP_PP_PATH is set!", color="OKGREEN")
    else:
        var = input("Please enter the VASP_PP_PATH environmen variable:\n")
        if len(var) == 0:
            cprint("Failed to set VASP_PP_PATH env. Abort.", color="FAIL")
            return False
        else:
            os.environ["VASP_PP_PATH"] = str(var)
            print(f"Set VASP_PP_PATH to {var}")
    return True


def demo_test():
    """Use low-level commands to execute VASP on a simple H2 structure
    check whether parsing of vasprun.xml, OUTCAR and vasp.out are supported
    """
    print("Running test example")
    atoms = molecule("H2", pbc=True, vacuum=4)
    with tempfile.TemporaryDirectory() as tmpdir:
        calc = VaspInteractive(
            nsw=0,
            istart=0,
            xc="pbe",
            directory=tmpdir,
            # directory="test1"
        )
        # Low level calculator interfacing

        with calc._txt_outstream() as out:
            try:
                # In some cases the vasp binary spin up can take longer
                with time_limit(30):
                    calc._run(atoms, out=out)
            except Exception as e:
                print("Running VASP encountered problem ", e)
                print(traceback.format_exc(limit=2))
                pass

        pid = calc.process.pid

        # Check vasprun.xml
        try:
            vrun = read(calc._indir("vasprun.xml"))
            # In newer version of ase the vasprun parsing is enhanced
            if vrun.calc is not None:
                vasprun_ok = True
            else:
                vasprun_ok = False
        except Exception as e:
            vasprun_ok = False

        # Check OUTCAR. In some cases the OUTCAR may be non-existing
        try:
            outcar_lines = open(calc._indir("OUTCAR"), "r").readlines()
        except Exception:
            outcar_lines = []
        outcar_ok = False
        cond = 0
        matches = [
            "FORCE on cell",
            "TOTAL-FORCE",
            "FREE ENERGIE OF THE ION-ELECTRON SYSTEM",
            "VOLUME and BASIS-vectors",
            "E-fermi",
        ]
        for i, line in enumerate(outcar_lines):
            if any([m in line for m in matches]):
                print(line)
                cond += 1
        if cond >= 5:
            outcar_ok = True

        # Check vaspout
        try:
            vaspout_lines = calc._txt_to_handler().readlines()
        except Exception:
            vaspout_lines = []
        # print(vaspout_lines)
        vaspout_ok = False
        for line in vaspout_lines:
            if "FORCES" in line:
                vaspout_ok = True
                break
        # Low level kill to prevent any issue with STOPCAR etc.
        subprocess.run(["kill", "-9", str(pid)])

    print("Single point calculation finished. Checking output file parsing.")
    status = None
    if vasprun_ok:
        cprint(f"vasprunxml: OK", color="OKGREEN")
    else:
        cprint(f"vasprunxml: FAIL", color="FAIL")
    if outcar_ok:
        cprint(f"OUTCAR: OK", color="OKGREEN")
    else:
        cprint(f"OUTCAR: FAIL", color="FAIL")
    if vaspout_ok:
        cprint(f"VASP raw output: OK", color="OKGREEN")
    else:
        cprint(f"VASP raw output: FAIL", color="FAIL")

    if not vaspout_ok:
        cprint(
            "VaspInteractive may not be compatible with your VASP setup. Please refer to the README for details.",
            color="FAIL",
        )
        status = "incompatible"
    elif not outcar_ok:
        cprint(
            (
                "VaspInteractive can read from raw output but OUTCAR is incomplete. "
                "Only the energy and forces are read in this case. "
                "Please refer to the README for details"
            ),
            color="WARNING",
        )
        status = "minimal support"
    elif not vasprun_ok:
        cprint(
            (
                "vasprun.xml is incomplete. "
                "VaspInteractive should still work but "
                "you're welcome to submit issues."
            ),
            color="OKBLUE",
        )
        status = "partial pass"
    else:
        cprint("All test pass! Enjoy coding.", color="OKGREEN")
        status = "all pass"

    # Output the total status
    print("#" * 80)
    print(f"Test result: {status}")
    print("#" * 80)


if __name__ == "__main__":
    if not vasp_env_test():
        sys.exit(1)
    demo_test()
