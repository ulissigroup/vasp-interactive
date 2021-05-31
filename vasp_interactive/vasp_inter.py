"""An I/O stream-based VASP calculator
Provides additional bug fix to
https://gitlab.com/ase/ase/-/blob/master/ase/calculators/vasp/interactive.py
May be merged with upstream by the end of day.
"""
from subprocess import Popen, PIPE
from contextlib import contextmanager

from ase.calculators.calculator import Calculator, ReadError
from ase.calculators.vasp import Vasp
from ase.io import read
from ase.io.vasp import write_vasp

from ase.calculators.vasp.vasp import check_atoms

import time
import os
import sys
import re
import numpy as np


class VaspInteractive(Vasp):
    """I/O stream-based VASP calculator.
    Can be used to speed up structural relaxation when using ASE optimizers.
    The atomic positions for each ionic step is passed to VASP via stdin.
    Currently do not support change of unit cell.
    """

    name = "VaspInteractive"
    # In principle the implemented properties can be anything
    # single point Vasp is capable of
    # Currently limits to e, F and S
    #     implemented_properties = ["energy", "forces", "stress"]
    implemented_properties = Vasp.implemented_properties
    mandatory_input = {
        #         "potim": 0.0,
        "ibrion": -1,
        #         "iwavpr": 11,
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
        txt="vasp.out",
        **kwargs,
    ):
        """Initialize the calculator object like the normal Vasp object.
        Additional attributes:
            `self.process`: Popen instance to run the VASP calculation
        At this stage, `VaspInteractive` does not allow restart
        """

        # Add the mandatory keywords
        for kw, val in VaspInteractive.mandatory_input.items():
            if kw in kwargs and val != kwargs[kw]:
                raise ValueError(
                    "Keyword {} cannot be overridden! "
                    "It must have have value {}, but {} "
                    "was provided instead.".format(kw, val, kwargs[kw])
                )
        kwargs.update(VaspInteractive.mandatory_input)

        for kw, val in VaspInteractive.default_input.items():
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
            **kwargs,
        )
        # VaspInteractive can take 1 Popen process to track the VASP job
        self.process = None

        # Make command a list of args for Popen
        cmd = self.make_command(self.command)
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._args = cmd

        # Ionic steps counter. Note this number will be 1 more than that in ase.optimize
        self.steps = 0
        # Is the relaxation finished?
        self.final = False

        return

    def reset(self):
        """Rewrite the parent reset function.
        The calculator is reset only if user set it to final
        """
        super(VaspInteractive, self).reset()
        if self.final:
            self._force_kill_process()
            self.steps = 0
            self.final = False

    def _ensure_directory(self):
        """Makesure self.directory exists, if not use `os.makedirs`"""
        # Create the folders where we write the files, if we aren't in the
        # current working directory.
        if self.directory != os.curdir and not os.path.isdir(
            self.directory
        ):
            os.makedirs(self.directory)

    @contextmanager
    def _txt_outstream(self):
        """Overwrites the parent Vasp._txt_outstream, so that the file io uses append mode
        Custom function for opening a text output stream. Uses self.txt
        to determine the output stream, and accepts a string or an open
        writable object.
        If a string is used, a new stream is opened, and automatically closes
        the new stream again when exiting.
        """

        txt = self.txt
        open_and_close = False  # Do we open the file?

        if txt is None:
            # Suppress stdout
            out = None
        else:
            if isinstance(txt, str):
                if txt == "-":
                    # subprocess.call redirects this to stdout
                    out = None
                else:
                    # Open the file in the work directory
                    self._ensure_directory()
                    txt = self._indir(txt)
                    # We wait with opening the file, until we are inside the
                    # try/finally
                    open_and_close = True
            elif hasattr(txt, "write"):
                out = txt
            else:
                raise RuntimeError(
                    "txt should either be a string"
                    "or an I/O stream, got {}".format(txt)
                )

        try:
            if open_and_close:
                # For interactive vasp mode, the io stream is in append mode
                # Using multiple stream context would not affect the performance
                out = open(txt, "a")
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
            raise RuntimeError(
                "VaspInteractive does not have the VASP process."
            )

    def _stdout(self, text, out=None):
        """ """
        if out is not None:
            out.write(text)

    def _run(self, atoms, out):
        """Overwrite the Vasp._run method
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
            #             print("Still running", self.process, self.process.poll())
            self._stdout("Inputting positions...\n", out=out)
            for atom in atoms.get_scaled_positions():
                self._stdin(
                    " ".join(map("{:19.16f}".format, atom)), out=out
                )

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


    def close(self):
        """Soft stop approach for the stream-based VASP process
        Works by writing STOPCAR file and runs two dummy scf cycles
        """
        if self.process is None:
            return

        #         print("Waiting to close ", self.process)
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
            #             print(
            #                 "Two consequetive runs of vasp for STOPCAR to work",
            #                 self.process,
            #                 self.process.poll(),
            #             )
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
            # No need to calculate, calculator has stored calculated results
            return

        # Currently VaspInteractive only handles change of positions (MD-like)
        if any([p in system_changes for p in ("numbers", "cell")]):
            # If not started yet,
            #             if self.steps == 0:
            #                 self.close()
            if self.process is not None:
                raise NotImplementedError(
                    (
                        "VaspInteractive does not support change of "
                        "chemical formula or lattice parameters. "
                        "Please create a new calculator instance or use standard Vasp calculator"
                    )
                )

        # TODO: add the out component
        with self._txt_outstream() as out:
            self._run(self.atoms, out=out)
            self.steps += 1
        
        # Use overwritten `read_results` method
        self.read_results()
        return
    
    def read_results(self):
        """Overwrites the `read_results` from parent class.
        In the interactive mode, after each ionic SCF cycle,
        only the OUTCAR content is written, while vasprun.xml
        is completed after user input. The results are read as 
        much as possible from the OUTCAR file.
        """
        # Temporarily load OUTCAR into memory
        outcar = self.load_file('OUTCAR')

        # vasprun.xml is only valid iteration when atoms finalized
        calc_xml = None
        xml_results = None
        if self.final:
            try:
                calc_xml = self._read_xml()
                xml_results = calc_xml.results
            except ReadError:
                # The xml file is not complete, try using OUTCAR only
                pass

        # Fix sorting
        if xml_results:
            xml_results['forces'] = xml_results['forces'][self.resort]
            self.results.update(xml_results)

        # OUTCAR handling part
        self.converged = self.read_convergence(lines=outcar)
        self.version = self.read_version()
        
        # Energy and magmom have multiple return values
        if "free_energy" not in self.results.keys():
            try:
                energy_free, energy_zero = self.read_energy(lines=outcar)
                self.results.update(dict(free_energy=energy_free,
                                         energy=energy_zero))
            except Exception:
                pass
            
        
        if "magmom" not in self.results.keys():
            try:
                magmom, magmoms = self.read_mag(lines=outcar)
                self.results.update(dict(magmom=magmom, 
                                         magmoms=magmoms))
            except Exception:
                pass
        
        # Missing properties that are name-consistent so can use dynamic function loading
        properties = ["forces", "stress", "fermi", "nbands", "dipole"]
        for prop in properties:
            if prop not in self.results.keys():
                try:
                    # use read_xxx method to parse outcar
                    result = getattr(self, f"read_{prop}")(lines=outcar)
                    self.results[prop] = result
                except Exception:
                    # Do not add the key
                    pass
                
        # Manunal old keywords handling
        self.spinpol = self.read_spinpol(lines=outcar)
#         self._set_old_keywords()

        # Store the parameters used for this calculation
        self._store_param_state()


    def read_all_iterations(self):
        """Parse the ionic & electronic scf cycles from OUTCAR files.
        Ideas taken from Vasp.read_number_of_iterations and Vasp.read_number_of_ionic_steps
           returns
           `n_ion_scf`: number of ionic steps
           `n_elec_scf`: list of electronic scf steps with length of n_ion
           
        """
        with self.load_file_iter("OUTCAR") as lines:
            n_ion_scf, n_elec_scf = parse_outcar_iterations(lines)
        return n_ion_scf, n_elec_scf
        
    def read_run_time(self):
        """Parse processing time from OUTCAR.
        returns (cpu_time, wall_time)
        If calculation is not finished, both are None
        """
        with self.load_file_iter("OUTCAR") as lines:
            cpu_time, wall_time = parse_outcar_time(lines)
        return cpu_time, wall_time
    
    def finalize(self):
        """Stop the stream calculator and finalize"""
        self._force_kill_process()
        self.final = True
        return
            
    def __enter__(self):
        """Reset everything upon entering the context"""
        self.reset()
        return self

    def __exit__(self, type, value, traceback):
        """Exiting the context manager and reset process"""
        self.finalize()
        return

    def _force_kill_process(self):
        """Try to kill the process by soft stop. If fails, force killing it"""
        try:
            self.close()
        # Runtime exceptions can occur when file id missing etc
        # Normally we don't want this behavior but just for backward-compatibility
        except Exception as e:
            # Do not use self.txt as output stream as it may not exist at this moment
            print(
                (
                    f"Trying to close the VASP stream but encountered error: \n"
                    f"{e}\n"
                    "Will now force closing the VASP process. "
                    "The OUTCAR and vasprun.xml outputs may be incomplete"
                ),
                file=sys.stderr,
            )
            if self.process is not None:
                if self.process.poll() is None:
                    self.process.kill()
        return

    def __del__(self):
        """Explicit deconstruction, kill the process with no mercy"""
        self._force_kill_process()
        return

    
# Following are functions parsing OUTCAR files which are not present in parent
# Vasp calculator but can he helpful for job diagnosis

def parse_outcar_iterations(lines):
    """Read the whole iteration information (ionic + electronic) from OUTCAR lines
    """
    n_ion_scf = 0
    n_elec_scf = []

    for line in lines:
        if '- Iteration' in line:
            ni_, ne_ = list(map(int, re.findall(r'\d+', line)))
            if ni_ > n_ion_scf:
                n_ion_scf = ni_
                n_elec_scf.append(ne_)
            else:
                n_elec_scf[ni_ - 1] = ne_
    n_elec_scf = np.array(n_elec_scf)
    return n_ion_scf, n_elec_scf

def parse_outcar_time(lines):
    """Parse the cpu and wall time from OUTCAR.
    The mismatch between wall time and cpu time represents 
    the turn-around time in VaspInteractive
    
    returns (cpu_time, wall_time)
    if the calculation is not finished, both will be None
    """
    cpu_time = None
    wall_time = None
    for line in lines:
        if "Total CPU time used (sec):" in line:
            cpu_time = float(line.split(":")[1].strip())
        if "Elapsed time (sec):" in line:
            wall_time = float(line.split(":")[1].strip())
    return cpu_time, wall_time