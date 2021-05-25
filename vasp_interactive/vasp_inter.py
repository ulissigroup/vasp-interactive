"""An I/O stream-based VASP calculator
   Provides additional bug fix to 
   https://gitlab.com/ase/ase/-/blob/master/ase/calculators/vasp/interactive.py
   May be merged with upstream by the end of day.
"""
from subprocess import Popen, PIPE

from ase.calculators.calculator import Calculator
from ase.calculators.vasp import Vasp
from ase.io import read

# from ase.calculators.vasp.create_input import GenerateVaspInput

import time
import os
import sys


class VaspInteractive(Vasp):
    name = "VaspInteractive"
    # In principle the implemented properties can be anything
    # single point Vasp is capable of
    # Currently limits to e, F and S
    implemented_properties = ["energy", "forces", "stress"]
    mandatory_input = {
        "potim": 0.0,
        "ibrion": -1,
        "iwavpr": 11,
        "interactive": True,
    }
    # Enforce the job to be relaxation
    default_input = {
        "nsw": 2000,
    }

    def __init__(
        self,
        atoms=None,
        directory=".",
        label="vasp-interactive",
        ignore_bad_restart_file=Calculator._deprecated,
        command=None,
        txt="vasp-interactive.out",
        print_log=False,
        **kwargs
    ):
        """Initialize the calculator object like the normal Vasp object.
        Additional attributes:
            `self.process`: Popen instance to run the VASP calculation
        At this stage, `VaspInteractive` does not allow restart
        """

        # Add the mandatory keywords
        for kw, val in self.mandatory_input.items():
            if kw in kwargs and val != kwargs[kw]:
                raise ValueError(
                    "Keyword {} cannot be overridden! "
                    "It must have have value {}, but {} "
                    "was provided instead.".format(kw, val, kwargs[kw])
                )
        kwargs.update(self.mandatory_input)

        for kw, val in self.default_input.items():
            if kw not in kwargs:
                kwargs[kw] = val

        super(VaspInteractive, self).__init__(
            atoms=atoms,
            directory=directory,
            label=label,
            ignore_bad_restart_file=ignore_bad_restart_file,
            command=command,
            txt=txt,
            restart=None,
            **kwargs
        )
        # VaspInteractive can take 1 Popen process to track the VASP job
        self.process = None
        self.print_log = print_log

        # The output context is slightly different from original Vasp calculator though
        if self.txt is None:
            # If txt is None, output stream will be supressed
            self.txt = None
        elif self._txt == "-":
            # If txt is '-' the output will be sent through stdout
            self.txt = sys.stdout
        elif isinstance(self.txt, str):
            # If txt is a string a file will be opened
            self.txt = open(self.txt, "a")
        else:
            # No change, but make sure self.txt has write
            if not hasattr(self.txt, "write"):
                raise AttributeError(
                    "If Vasp.txt should be either a string, PurPath or IO stream."
                )
            pass

        # Make command a list of args for Popen
        cmd = self.make_command(self.command)
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._args = cmd

        return

    def _stdin(self, text, ending="\n"):
        """Pass input to the stdin of VASP process
        single line
        """
        if self.txt is not None:
            self.txt.write(text + ending)
        if self.print_log:
            print(text, end=ending)
        if self.process is not None:
            self.process.stdin.write(text + ending)
            # ASE only supports py3 now, no need for py2 compatibility
            self.process.stdin.flush()
        else:
            raise RuntimeError("VaspInteractive does not have the VASP process.")

    def _stdout(self, text):
        if self.txt is not None:
            self.txt.write(text)
        if self.print_log:
            print(text, end="")

    def _run_vasp(self, atoms):
        if self.process is None:
            print("Initialize")
            stopcar = os.path.join(self.directory, "STOPCAR")
            if os.path.isfile(stopcar):
                os.remove(stopcar)
            self._stdout("Writing VASP input files\n")
            self.initialize(atoms)
            self.write_input(atoms)
            self._stdout("Starting VASP for initial step...\n")
            # Drop py2 support
            self.process = Popen(
                        self._args,
                        stdout=PIPE,
                        stdin=PIPE,
                        stderr=PIPE,
                        cwd=self.directory,
                        universal_newlines=True,
                        bufsize=0,
                    )
        else:
            print("Still running", self.process, self.process.poll())
            self._stdout("Inputting positions...\n")
            for atom in atoms.get_scaled_positions():
                self._stdin(" ".join(map("{:19.16f}".format, atom)))

        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                return

        # Extra condition, vasp exited with 0 (completed)
        # only happens at second call to _run_vasp when STOPCAR is present
        if self.process.poll() == 0:
            self._stdout("VASP successfully terminated")
            return
        else:
            # If we've reached this point, then VASP has exited without asking for
            # new positions, meaning it either exited without error unexpectedly,
            # or it exited with an error. Either way, we need to raise an error.

            raise RuntimeError(
                "VASP exited unexpectedly with exit code {}"
                "".format(self.process.poll())
            )

    def close(self):
        if self.process is None:
            return

        print("Waiting to close ", self.process)
        self._stdout("Attemping to close VASP cleanly\n")
        with open(os.path.join(self.directory, "STOPCAR"), "w") as stopcar:
            stopcar.write("LABORT = .TRUE.")

        # Following two calls to _run_vasp: 1 is to let vasp see STOPCAR and do 1 SCF
        # second is to exit the program
        self._run_vasp(self.atoms)
        self._run_vasp(self.atoms)
        # if runs to this stage, process.poll() should be 0
        print(
            "Two consequetive runs of vasp for STOPCAR to work",
            self.process,
            self.process.poll(),
        )
        while self.process.poll() is None:
            time.sleep(1)
        self._stdout("VASP has been closed\n")
        self.process = None
        return

    def calculate(
        self,
        atoms=None,
        properties=["energy"],
        system_changes=["positions", "numbers", "cell"],
    ):

        atoms = atoms.copy()
        self.results.clear()
        Calculator.calculate(self, atoms, properties, system_changes)

        if not system_changes:
            return

        if "numbers" in system_changes:
            self.close()

        self._run_vasp(atoms)

        print("Before reading OUTCAR")
        max_retry = 1
        new = None
        for i in range(max_retry):
            try:
                new = read(
                    os.path.join(self.directory, "OUTCAR"),
                )
                # index=-1)
            except Exception as e:
                print(e)
                # Just continue the look
                print("Failed OUTCAR read for time {}".format(i))
                time.sleep(0.5)
                outcar = os.path.join(self.directory, "OUTCAR")
                # print(read(outcar))
                with open(outcar, "r") as fd:
                    print("OUTCAR looks like this (last 25 lines)")
                    lines = fd.readlines()
                    print("".join(lines[-25:]))
                    # for l in lines: print(l)
                    content = "".join(lines)
                    # from ase.io.vasp import read_vasp_out
                    # from io import StringIO
                    # from tempfile import NamedTemporaryFile
                    # print(read_vasp_out(outcar))
                    # from pymatgen.io.vasp import Outcar
                    # ot = Outcar(outcar)
                    # print(ot)

                with open("OUTCAR.error", "w") as fd:
                    print("Writing to OUTCAR.error")
                    fd.write("".join(lines))
                # tmp = NamedTemporaryFile()
                # tmp.write(content)
                # print(read_vasp_out("OUTCAR.error"))
                pass
        #         new = read(os.path.join(self.directory, 'vasprun.xml'), index=-1)
        # print("Finish reading outcar")
        if new:
            self.results = {
                "free_energy": new.get_potential_energy(force_consistent=True),
                "energy": new.get_potential_energy(),
                "forces": new.get_forces()[self.resort],
            }
        else:
            from pymatgen.io.vasp import Outcar
            import numpy as np

            # Dirty patch, not using os.path
            # from https://gist.github.com/gVallverdu/0e232988f32109b5dc6202cf193a49fb
            ot = Outcar("OUTCAR")
            forces = ot.read_table_pattern(
                header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
                row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
                footer_pattern=r"\s--+",
                postprocess=lambda x: float(x),
                last_one_only=False,
            )
            forces = np.array(forces[-1])
            self.results = {
                "free_energy": ot.final_energy,
                "energy": ot.final_energy,
                "forces": forces[self.resort],
            }
        # print(self.resort)
        print(self.results)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


#     def __del__(self):
#         self.close()

# if __name__ == "__main__":
#     from ase import Atoms
#     from ase.build import molecule
#     from ase.optimize import BFGS, QuasiNewton
#     from ase.calculators.emt import EMT
#     import numpy as np
#     import os

#     os.environ["VASP_COMMAND"] = "mpirun -np 8 /opt/vasp.6.1.2_pgi_mkl_beef/bin/vasp_std"

#     # Delete all *CAR files to reproduce error
#     for f in next(os.walk("./"))[-1]:
#         if "CAR" in f:
#             os.unlink(f)

#     parent_calc = VaspInteractive(
#         istart=0,
#         algo="Fast",
#         prec="Normal",
#         isif=4,
#         ismear=0,
#         ispin=1,
#         ediff=1e-4,
#         xc="rpbe",
#         encut=300,
#         #     nsw=0,  # important that this is set higher than the number of times it will be called
#         kpts=(1, 1, 1),  # not important, just keeps it faster
#     )

#     # atoms = molecule("C6H6", vacuum=5)
#     # atoms.rattle(0.05)
#     # atoms.set_calculator(parent_calc)
#     with parent_calc:
#         d = 0.9575
#         h2 = Atoms(
#             "H2", positions=[(d, 0, 0), (0, 0, 0)], cell=[8, 8, 8], calculator=parent_calc
#         )

#         # dyn = BFGS(atoms)
#         dyn = BFGS(h2)
#         dyn.run(fmax=0.05)
