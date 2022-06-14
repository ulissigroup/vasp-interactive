"""An I/O stream-based VASP calculator
Provides additional bug fix to
https://gitlab.com/ase/ase/-/blob/master/ase/calculators/vasp/interactive.py
May be merged with upstream by the end of day.
"""
from subprocess import Popen, PIPE
from contextlib import contextmanager

from ase.calculators.calculator import Calculator, ReadError, CalculatorSetupError
from ase.calculators.vasp import Vasp
from ase.io import read
from ase.io.vasp import write_vasp

from ase.calculators.vasp.vasp import check_atoms
from warnings import warn

import time
import os
import sys
import psutil
import signal
import re
import numpy as np

DEFAULT_KILL_TIMEOUT = 60


class TimeoutException(Exception):
    """Simple class for timeout"""

    pass


@contextmanager
def time_limit(seconds):
    """Usage:
    try:
        with time_limit(60):
            do_something()
    except TimeoutException:
        raise
    """

    def signal_handler(signum, frame):
        raise TimeoutException("Timed out closing VASP process.")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def _find_mpi_process(pid, mpi_program="mpirun", vasp_program="vasp_std"):
    """Recursively search children processes with PID=pid and return the one
    that mpirun (or synonyms) are the main command.
    """
    allowed_names = set(["mpirun", "mpiexec", "orterun", "oshrun", "shmemrun"])
    allowed_vasp_names = set(["vasp_std", "vasp_gam", "vasp_ncl"])
    if mpi_program:
        allowed_names.add(mpi_program)
    if vasp_program:
        allowed_vasp_names.add(vasp_program)
    process_list = [psutil.Process(pid)]
    process_list.extend(process_list[0].children(recursive=True))
    mpi_candidates = []
    for proc in process_list:
        # print(proc, proc.name())
        if proc.name() in allowed_names:
            # is the mpi process's direct children are vasp_std?
            children = proc.children()
            if len(children) > 0:
                if children[0].name() in allowed_vasp_names:
                    mpi_candidates.append(proc)
    if len(mpi_candidates) > 1:
        warn(
            "More than 1 mpi processes are created. This may be a bug. I'll use the last one"
        )
    mpi_proc = mpi_candidates[-1]
    return mpi_proc


class VaspInteractive(Vasp):
    """I/O stream-based VASP calculator.
    Can be used to speed up structural relaxation when using ASE optimizers.
    The atomic positions for each ionic step is passed to VASP via stdin.
    Currently do not support change of unit cell.
    """

    name = "VaspInteractive"
    implemented_properties = Vasp.implemented_properties
    mandatory_input = {
        "potim": 0.0,
        "ibrion": -1,
        "iwavpr": 11,
        "interactive": True,
        # Disable stopping criteria but just rely on nsw
        "ediffg": 0,
    }
    # Enforce the job to run infinitely untill killed
    default_input = {
        "nsw": 2000,
        "isym": 0,
    }

    def __init__(
        self,
        atoms=None,
        directory=".",
        label="vasp-interactive",
        ignore_bad_restart_file=Calculator._deprecated,
        command=None,
        txt="vasp.out",
        allow_restart_process=True,
        allow_mpi_pause=True,
        allow_default_param_overwrite=True,
        cell_tolerance=1e-8,
        kill_timeout=DEFAULT_KILL_TIMEOUT,
        **kwargs,
    ):
        """Initialize the calculator object like the normal Vasp object.
        Additional attributes:
            `self.process`: Popen instance to run the VASP calculation
            `allow_restart_process`: if True, will restart the VASP process if it exits before user program ends
            `allow_mpi_pause`: If disabled, do not interfere with the vasp program but let system load balancing handle CPU requests.
            `allow_default_param_overwrite`: If True, use mandatory input parameters to overwrite (but give warnings)
        """

        # Add the mandatory keywords
        for kw, val in VaspInteractive.mandatory_input.items():
            if kw in kwargs and val != kwargs[kw]:
                if allow_default_param_overwrite:
                    warn(
                        f"You have provided {kw} with value {kwargs[kw]}. "
                        f"For VaspInteractive to run properly it needs to be {val}. "
                        "I will overwrite."
                    )
                else:
                    raise ValueError(
                        f"Keyword {kw} cannot be overridden! "
                        f"It must have have value {val}, but {kwargs[kw]} "
                        "was provided instead."
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
        self.allow_restart_process = allow_restart_process

        # Ionic steps counter. Note this number will be 1 more than that in ase.optimize
        self.steps = 0
        # Is the relaxation finished?
        self.final = False

        # If nsw=0 or 1, the user are using VaspInteractive as normal Vasp single point
        # can generate warning
        incar_nsw = self.int_params["nsw"]
        if incar_nsw in (0, 1):
            warn(
                (
                    f"You have set NSW={incar_nsw} in INCAR. "
                    "VaspInteractive will run as a normal single point calculator. "
                    "If this is what you want, ignore this warning."
                )
            )

        # Recommend to use isym=0 like in MD calculations, otherwise energy can be wrong
        incar_isym = self.int_params["isym"]
        if incar_isym > 0:
            # If user defines isym to another value, use at their risk
            warn(
                "It is recommended to use ISYM=0 to "
                "overcome convergence issue when input symmetry changes. \n"
                f"However, you provided ISYM={incar_isym}. "
                "In some cases the energy and forces can be wrong. "
                "Use such settings at your own risk."
            )

        # Cell tolerance parameter
        self.cell_tolerance = abs(cell_tolerance)
        if self.cell_tolerance > 1e-3:
            warn(
                f"Your cell tolerace of {self.cell_tolerance} is probably too high. "
                "Make sure your results make sense"
            )

        # Add pause function
        self.pause_mpi = allow_mpi_pause
        self.kill_timeout = kill_timeout
        return

    @property
    def incar_nsw(self):
        nsw_ = self.int_params["nsw"]
        if nsw_ == 0:
            nsw_ = 1
        return nsw_

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
                    out = sys.stdout
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
            raise RuntimeError("VaspInteractive does not have the VASP process.")

    def _stdout(self, text, out=None):
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
        # VASP process exited (possibly due to NSW limit reached), release the process handle
        # a bit messy conditions here but works
        if self.process is not None:
            if (self.process.poll() == 0) and self.allow_restart_process:
                pid = self.process.pid
                self.process = None
                self._stdout(
                    "It seems your VASP process exited normally. I'll retart a new one.",
                    out=out,
                )
                warn(
                    (
                        f"VASP process (pid={pid}) exits normally but new positions are still provided. "
                        "A new VASP process will be started. "
                        "To supress this warning, you may want to increase the NSW number in your settings."
                    )
                )

        if self.process is None:
            # Delete STOPCAR left by an unsuccessful run
            stopcar = self._indir("STOPCAR")

            if os.path.isfile(stopcar):
                os.remove(stopcar)
            self._stdout("Writing VASP input files\n", out=out)
            self.initialize(atoms)
            self.write_input(atoms)
            self._stdout("Starting VASP for initial step...\n", out=out)
            # Dynamic generation of command args
            command = self.make_command(self.command)
            self.process = Popen(
                command,
                shell=True,
                stdout=PIPE,
                stdin=PIPE,
                stderr=PIPE,
                cwd=self.directory,
                universal_newlines=True,
                bufsize=0,
            )
            self.steps = 0
        else:
            # Whenever at this point, VASP interactive asks the input
            # write the current atoms positions to the stdin
            retcode = self.process.poll()
            if retcode is None:
                self._stdout("Inputting positions...\n", out=out)
                # ase atoms --> self.sort --> write to position for VASP
                # use scaled positions wrap back to cell
                for atom in atoms.get_scaled_positions()[self.sort]:
                    self._stdin(" ".join(map("{:19.16f}".format, atom)), out=out)
            else:
                # The vasp process stops prematurely
                raise RuntimeError(
                    (
                        f"The VASP process has exited with code {retcode} but "
                        "you're still providing new atom positions. "
                        "If the return code is 0, you can try to set allow_restart_process=True "
                        "to enable auto restart of the VASP process."
                    )
                )

        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text, out=out)
            # Read vasp version from stdio, warn user of VASP5 issue
            if self._read_vasp_version_stream(text):
                if _int_version(self.version) < 6:
                    warn(
                        (
                            "Some builds of VASP 5.x may have issue generating output blocks. "
                            "See this issue for details https://github.com/ulissigroup/vasp-interactive/issues/6. "
                            "If you encounter similar error messages, try using VASP version > 6 if available."
                        )
                    )
            if "POSITIONS: reading from stdin" in text:
                return

        # Extra condition, vasp exited with 0 (completed)
        # can be 2 situations: completed job due to STOPCAR
        # or nsw is reached
        if self.process.poll() == 0:
            self._stdout("VASP terminated normally\n", out=out)
            # Update Aug. 13 2021
            # Since we explicitly added the check for self.steps in self.calculate
            # the following scenario should not happen
            if self.steps > self.incar_nsw:
                self._stdout(
                    (
                        "However the maximum ionic iterations have been reached. "
                        "Consider increasing NSW number in your calculation."
                    ),
                    out=out,
                )
                raise RuntimeError(
                    (
                        "VASP process terminated normally but "
                        "your current ionic steps exceeds maximum allowed number. "
                        "Consider increase your NSW value in calculator setup, "
                        "or set allow_process_restart=True"
                    )
                )

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
        #TODO: add a time-out option
        """
        # Make sure the MPI process is awaken before the termination
        if self.pause_mpi:
            self._resume_calc()

        if self.process is None:
            return
        elif self.process.poll() is not None:
            # For whatever reason the vasp process stops prematurely (possibly too small nsw)
            # do a clean up
            retcode = self.process.poll()
            with self._txt_outstream() as out:
                self._stdout(f"VASP exited with code {retcode}.", out=out)
            self.process = None
            return
        else:
            with self._txt_outstream() as out:
                self._stdout("Attemping to close VASP cleanly\n", out=out)
                stopcar = self._indir("STOPCAR")
                with open(stopcar, "w") as fd:
                    fd.write("LABORT = .TRUE.")

                # Following two calls to _run_vasp: 1 is to let vasp see STOPCAR and do 1 SCF
                # second is to exit the program
                # Program may have ended before we can write, if so then cancel the stdin
                for i in range(2):
                    self._run(self.atoms, out=out)
                    if self.process.poll() is not None:
                        self._stdout(
                            f"VASP exited with code {self.process.poll()}.", out=out
                        )
                        self.process = None
                        return
                # TODO: the endless waiting cycle is hand-waving
                # consider add a timeout function
                while self.process.poll() is None:
                    time.sleep(1)
                self._stdout("VASP has been closed\n", out=out)
                self.process = None
            return

    def _pause_calc(self, sig=signal.SIGTSTP):
        """Pause the vasp processes by sending SIGTSTP to the master mpirun process"""
        if not self.process:
            return
        pid = self.process.pid
        mpi_process = _find_mpi_process(pid)
        if mpi_process is None:
            warn("Cannot find the mpi process. Will not send stop signal to mpi.")
            return
        mpi_process.send_signal(sig)
        return

    def _resume_calc(self, sig=signal.SIGCONT):
        """Resume the vasp processes by sending SIGCONT to the master mpirun process"""
        if not self.process:
            return
        pid = self.process.pid
        mpi_process = _find_mpi_process(pid)
        if mpi_process is None:
            warn("Cannot find the mpi process. Will not send continue signal to mpi.")
            return
        mpi_process.send_signal(sig)
        return

    @contextmanager
    def pause(self):
        """Context wrapper for pause / resume MPI process. To avoid accidentally forgetting the resume.
        Usage:
        with calc.pause():
            # Do some CPU expensive job here
        This context wrapper can be disabled by setting self.pause_mpi = False
        """
        if self.pause_mpi:
            try:
                self._pause_calc()
                yield
            finally:
                self._resume_calc()
        else:
            # Do nothing
            yield

    @contextmanager
    def _ensure_mpi(self):
        """Opposite to the pause context manager,
        _ensure_mpi ensures the MPI process is running only within the block
        """
        try:
            self._resume_calc()
            yield
        finally:
            self._pause_calc()

    def check_state(self, atoms, tol=1e-15):
        """Modified check_state method to allow separate check for cell tolerance"""
        old_system_changes = super(VaspInteractive, self).check_state(atoms, tol=1e-15)
        if ("cell" in old_system_changes) and (self.atoms is not None):
            max_cell_change = np.max(np.abs(atoms.cell - self.atoms.cell))
            # Do not set cell change
            if max_cell_change < self.cell_tolerance:
                old_system_changes = [sc for sc in old_system_changes if sc != "cell"]
        return old_system_changes

    def calculate(
        self,
        atoms=None,
        properties=["energy"],
        system_changes=["positions", "numbers", "cell"],
    ):
        check_atoms(atoms)

        if hasattr(self, "system_changes") and self.system_changes is not None:
            system_changes = self.system_changes
            self.system_changes = None

        if not system_changes:
            # No need to calculate, calculator has stored calculated results
            return

        # Currently VaspInteractive only handles change of positions (MD-like)
        if "numbers" in system_changes:
            if self.process is not None:
                raise NotImplementedError(
                    (
                        "VaspInteractive does not support change of chemical formula. "
                        "Please create a new calculator instance or use standard Vasp calculator"
                    )
                )
        elif "cell" in system_changes:
            if self.process is not None:
                raise NotImplementedError(
                    (
                        "VaspInteractive does not support change of lattice parameters. "
                        "Set VaspInteractive.cell_tolerance to a higher value if you think it's caused by round-off error. "
                        "Otherwise, please create a new calculator instance or use standard Vasp calculator"
                    )
                )

        self.clear_results()
        if atoms is not None:
            self.atoms = atoms.copy()

        with self._txt_outstream() as out:
            self._run(self.atoms, out=out)
            self.steps += 1
            # special condition: job runs with nsw limit reached.
            # In the interactive mode, VASP won't exit until steps > nsw
            # this can be problematic if the next user input is a different position
            # so we simply run a dummy step using current positions to gracefully terminate VASP
            if self.steps >= self.incar_nsw:
                self._run(self.atoms, out=out)

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
        outcar = self.load_file("OUTCAR")

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
            xml_results["forces"] = xml_results["forces"][self.resort]
            self.results.update(xml_results)

        # OUTCAR handling part
        self.converged = self.read_convergence(lines=outcar)
        self.version = self.read_version()
        if self.version is not None:
            if _int_version(self.version) < 6:
                warn(
                    (
                        "Some builds of VASP 5.x may have issue generating output blocks. "
                        "See this issue for details https://github.com/ulissigroup/vasp-interactive/issues/6. "
                        "If you encounter similar error messages, try using VASP version > 6 if available."
                    )
                )

        # Energy and magmom have multiple return values
        if "free_energy" not in self.results.keys():
            try:
                energy_free, energy_zero = self.read_energy(lines=outcar)
                self.results.update(dict(free_energy=energy_free, energy=energy_zero))
            except Exception:
                pass

        if "magmom" not in self.results.keys():
            try:
                magmom, magmoms = self.read_mag(lines=outcar)
                self.results.update(dict(magmom=magmom, magmoms=magmoms))
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

    def _read_vasp_version_stream(self, line):
        """Read vasp version from streamed output lines"""
        if " vasp." in line:
            self.version = line[len(" vasp.") :].split()[0]
            return True
        else:
            return False

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
            with time_limit(self.kill_timeout):
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

    def __copy__(self):
        """Overwrite the shallow copy method
        this will be the same behavior as standard shallow copy,
        while we manually put self.process to None.

        Note shallow copying the calculator won't copy the list / dict attributes
        """
        from warnings import warn

        new = type(self)()
        new.__dict__.update(self.__dict__)
        new.process = None
        warn(
            (
                "Due to thread safety, copy of the VaspInteractive calculator "
                "does not contain the reference to the VASP process."
            )
        )
        return new

    def __deepcopy__(self, memo):
        """Overwrite the deepcopy method, use deepcopy to copy every attribute except process
        see https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object/40484215"""
        from warnings import warn
        from copy import deepcopy

        cls = self.__class__
        new = cls.__new__(cls)
        memo[id(self)] = new
        for k, v in self.__dict__.items():
            if k != "process":
                setattr(new, k, deepcopy(v, memo))
            else:
                setattr(new, k, None)
        warn(
            (
                "Due to thread safety, deepcopy of the VaspInteractive calculator "
                "does not contain the reference to the VASP process"
            )
        )
        return new


# Following are functions parsing OUTCAR files which are not present in parent
# Vasp calculator but can he helpful for job diagnosis


def parse_outcar_iterations(lines):
    """Read the whole iteration information (ionic + electronic) from OUTCAR lines"""
    n_ion_scf = 0
    n_elec_scf = []

    for line in lines:
        if "- Iteration" in line:
            ni_, ne_ = list(map(int, re.findall(r"\d+", line)))
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


def _int_version(version_string):
    """Get int version string"""
    major = int(version_string.split(".")[0])
    return major
