"""An I/O stream-based VASP calculator
   Provides additional bug fix to 
   https://gitlab.com/ase/ase/-/blob/master/ase/calculators/vasp/interactive.py
   May be merged with upstream by the end of day.
"""
from subprocess import Popen, PIPE
from contextlib import contextmanager

from ase.calculators.calculator import Calculator
from ase.calculators.vasp import Vasp
from ase.io import read

from ase.calculators.vasp.vasp import check_atoms


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

        # Make command a list of args for Popen
        cmd = self.make_command(self.command)
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._args = cmd

        return
    
    def _ensure_directory(self):
        """Makesure self.directory exists, if not use `os.makedirs`
        """
        # Create the folders where we write the files, if we aren't in the
        # current working directory.
        if self.directory != os.curdir and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    
    @contextmanager
    def _txt_outstream(self):
        """Overwrites the parent Vasp._txt_outstream, so that the file io uses append mode
        Custom function for opening a text output stream. Uses self.txt
        to determine the output stream, and accepts a string or an open
        writable object.
        If a string is used, a new stream is opened, and automatically closes
        the new stream again when exiting.

        Examples:
        # Pass a string
        calc.txt = 'vasp.out'
        with calc.txt_outstream() as out:
            calc.run(out=out)   # Redirects the stdout to 'vasp-interactive.out'

        # Use an existing stream
        mystream = open('vasp.out', 'w')
        calc.txt = mystream
        with calc.txt_outstream() as out:
            calc.run(out=out)
        mystream.close()

        # Print to stdout
        calc.txt = '-'
        with calc.txt_outstream() as out:
            calc.run(out=out)   # output is written to stdout
        """

        txt = self.txt
        open_and_close = False  # Do we open the file?

        if txt is None:
            # Suppress stdout
            out = None
        else:
            if isinstance(txt, str):
                if txt == '-':
                    # subprocess.call redirects this to stdout
                    out = None
                else:
                    # Open the file in the work directory
                    self._ensure_directory()
                    txt = self._indir(txt)
                    # We wait with opening the file, until we are inside the
                    # try/finally
                    open_and_close = True
            elif hasattr(txt, 'write'):
                out = txt
            else:
                raise RuntimeError('txt should either be a string'
                                   'or an I/O stream, got {}'.format(txt))

        try:
            if open_and_close:
                # For interactive vasp mode, the io stream is in append mode
                # Using multiple stream context would not affect the performance
                out = open(txt, 'a')
            yield out
        finally:
            if open_and_close:
                out.close()


    def _stdin(self, text, out=None, ending="\n"):
        """Write single line text to `self.process.stdin`
           if `out` provided, the same input is write to `out` as well.
           `text` should be a single line input without ending char
        """
        if out is not None:
            out.write(text + ending)
        if self.process is not None:
            self.process.stdin.write(text + ending)
            # ASE only supports py3 now, no need for py2 compatibility
            self.process.stdin.flush()
        else:
            raise RuntimeError("VaspInteractive does not have the VASP process.")

    def _stdout(self, text, out=None):
        """ 
        """
        if out is not None:
            out.write(text)

    def _run(self, atoms, out):
        """ Overwrite the Vasp._run method
            Running pipe-based vasp job
            `out` is the io stream determined by `_txt_outstream()`
            Logic: 
            - If the vasp process not present (either not started or restarted):
                make inputs and run until stdout captures request for new pos input
            - Else the vasp process has started and asks for input (2+ ionic step)
                write new positions to stdin
            - If the process continues without asking input, check the process return code
              if returncode != 0 then there is an error
        """
        if self.process is None:
            # Delete STOPCAR left by an unsuccessful run
            stopcar = self._indir("STOPCAR")
            if os.path.isfile(stopcar):
                os.remove(stopcar)
            self._stdout("Writing VASP input files\n", out=out)
            self.initialize(atoms)
            self.write_input(atoms)
            self._stdout("Starting VASP for initial step...\n", out=out)
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
            # Whenever at this point, VASP interactive asks the input
            # write the current atoms positions to the stdin
            print("Still running", self.process, self.process.poll())
            self._stdout("Inputting positions...\n", out=out)
            for atom in atoms.get_scaled_positions():
                self._stdin(" ".join(map("{:19.16f}".format, atom)), out=out)

        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text, out=out)
            if "POSITIONS: reading from stdin" in text:
                return

        # Extra condition, vasp exited with 0 (completed)
        # only happens at second call to _run_vasp when STOPCAR is present
        if self.process.poll() == 0:
            self._stdout("VASP successfully terminated\n", out=out)
            return
        else:
            # If we've reached this point, then VASP has exited without asking for
            # new positions, meaning it either exited without error unexpectedly,
            # or it exited with an error. Either way, we need to raise an error.

            raise RuntimeError(
                "VASP exited unexpectedly with exit code {}"
                "".format(self.process.poll())
            )

#     def _close_io(self):
#         """ Explicitly close io stream of self.txt
            
#         """
#         if hasattr(self.txt, "write"):
#             if self.txt.closed is False:
#                 print("closing io stream", self.txt)
#                 self.txt.close()
#         return
        
    def close(self):
        """Necessary step to stop the stream-based VASP process
           Works by writing STOPCAR file and runs two dummy scf cycles
        """
        if self.process is None:
            return

        print("Waiting to close ", self.process)
        with self._txt_outstream() as out:
            self._stdout("Attemping to close VASP cleanly\n", out=out)
            stopcar = self._indir("STOPCAR")
            with open(stopcar, "w") as fd:
                fd.write("LABORT = .TRUE.")

            # Following two calls to _run_vasp: 1 is to let vasp see STOPCAR and do 1 SCF
            # second is to exit the program
            self._run(self.atoms, out=out)
            self._run(self.atoms, out=out)
            # if runs to this stage, process.poll() should be 0
            print(
                "Two consequetive runs of vasp for STOPCAR to work",
                self.process,
                self.process.poll(),
            )
            # TODO: the endless waiting cycle is hand-waving
            # consider add a timeout function
            while self.process.poll() is None:
                time.sleep(1)
            self._stdout("VASP has been closed\n", out=out)
            self.process = None
        return

    def calculate(
        self,
        atoms=None,
        properties=["energy"],
        system_changes=["positions", "numbers", "cell"],
    ):
        # TODO: use base method to handle directory
        check_atoms(atoms)
        self.clear_results()
        if atoms is not None:
            self.atoms = atoms.copy()

        if not system_changes:
            return
        
#         # TODO: Doubt about this part
#         # maybe something like self.restart() ?
        if "numbers" in system_changes:
            self.close()
        
        # TODO: add the out component
        with self._txt_outstream() as out:
            self._run(self.atoms, out=out)

        print("Before reading OUTCAR")
        max_retry = 1
        new = None
        outcar = self._indir("OUTCAR")
        for i in range(max_retry):
            try:
                new = read(outcar, index=-1)
            except Exception as e:
                print(e)
                # Just continue the look
                print("Failed OUTCAR read for time {}".format(i))
                time.sleep(0.5)
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

                with open(self._indir("OUTCAR.error"), "w") as fd:
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
            # Dirty patch, not using os.path
            # from https://gist.github.com/gVallverdu/0e232988f32109b5dc6202cf193a49fb
            from pymatgen.io.vasp import Outcar
            import numpy as np

            ot = Outcar(self._indir("OUTCAR"))
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
        """Close the process and file operators
        """
        self.close()
#         self._close_io()

    def __del__(self):
        """Close the process and file operators
        """
        self.close()
#         self._close_io()