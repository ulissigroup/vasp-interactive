"""An I/O stream-based VASP calculator
Provides additional bug fix to
https://gitlab.com/ase/ase/-/blob/master/ase/calculators/vasp/interactive.py
May be merged with upstream by the end of day.
"""
import time
import signal
import traceback
import pickle
import os
import re
import psutil
import sys
from pathlib import Path
from warnings import warn
from subprocess import Popen, PIPE
from contextlib import contextmanager

import numpy as np
from ase.calculators.calculator import (
    Calculator,
    ReadError,
    CalculatorSetupError,
    all_properties,
    all_changes,
)
from ase.calculators.vasp.vasp import Vasp, check_atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator, SinglePointCalculator
from ase.calculators.socketio import SocketClient


from .utils import (
    DEFAULT_KILL_TIMEOUT,
    _int_version,
    time_limit,
    _find_mpi_process,
    _slurm_signal,
    _preprocess_mlff_outcar,
)
from .parse import (
    parse_outcar_iterations,
    parse_outcar_time,
    parse_vaspout_energy,
    parse_vaspout_forces,
)


class VPISocketClient(SocketClient):
    """Minor patch to SocketClient to handle parent calculator termination"""

    def attach_parent_calc(self, calc):
        self.parent_calc = calc
        return

    def close(self):
        if not self.closed:
            self.log("Close SocketClient")
            self.protocol.socket.close()
            # If there is a calculator associated, use finalize to terminate
            if hasattr(self, "parent_calc"):
                if (self.parent_calc is not None) and (
                    hasattr(self.parent_calc, "finalize")
                ):
                    if not self.parent_calc.final:
                        self.parent_calc.finalize()
            self.closed = True


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
        # Set ISIF tag to cature stress-induced energy change
        "isif": 3,
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
        command=None,
        txt="vasp.out",
        allow_restart_process=True,
        allow_mpi_pause=True,
        allow_default_param_overwrite=True,
        cell_tolerance=1e-8,
        kill_timeout=DEFAULT_KILL_TIMEOUT,
        parse_vaspout=True,
        use_socket=False,
        host="localhost",
        port=None,
        unixsocket=None,
        timeout=None,
        log=None,
        **kwargs,
    ):
        """Interactive mode calculator for VASP.

        `VaspInteractive` handles all compatible parameters with `ase.calculators.vasp.Vasp`.
        There are also additional parameters controlling the interactive mode and socket-I/O mode
        behavior:

        Parameters for interactive mode:
            `allow_restart_process`: if True, will restart the VASP process if it exits before user program ends
            `allow_mpi_pause`: whether the calculator can be paused.
            `allow_default_param_overwrite`: if True, use mandatory input parameters to overwrite (but give warnings)
            `cell_tolerance`: tolerance threshold for detecting cell change.
            `kill_timeout`: timeout in seconds before forcibly kill the vasp process
            `parse_vaspout`: Whether to parse vasp.out for incorrect energy and forces fields. Only relevant if using VASP 5.x

        Parameters for socket-I/O mode:
            `use_socket`: if True, attach a socket client to the calculator and self.run() method
                          will be available. Note you don't need to set use_socket when passing through
                          SocketIOCalculator
            `host`: hostname of the socket server running iPI protocol
            `port`: server port to connect to
            `unixsocket`: name of local unix socket
            `timeout`: socket I/O timeout
            `log`: logfile for the socket client. If None, write to stdout

        Different behaviors compared with original `ase.calculators.vasp.Vasp`:
        - VaspInteractive can be run using the context manager (the `with`-clause)
        - Copies of the calculator do not duplicate the underlying VASP process
        - self.command is by default set to calling the socket mode command line.
          The actual VASP command passed in __init__ is stored in self._vasp_command
        """

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
            command=command,
            txt=txt,
            restart=None,
            **kwargs,
        )
        # VaspInteractive can take 1 Popen process to track the VASP job
        # self.pid tracks if pid changes (useful for stopping slurm jobs)
        self.process = None
        self.pid = None
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
        # Tracker of mpi process (root proc of mpirun or slurm stepid etc)
        # after calling _find_mpi_process, this property will become a dict cont
        if self.pause_mpi:
            self.mpi_match = None
            self.mpi_state = None
        self.kill_timeout = kill_timeout
        self.parse_vaspout = parse_vaspout

        # VASP status flags
        self._xml_complete = None
        self._outcar_complete = None

        input_params = dict(
            label=label,
            command=command,
            txt=txt,
            allow_restart_process=allow_restart_process,
            allow_mpi_pause=allow_mpi_pause,
            allow_default_param_overwrite=allow_default_param_overwrite,
            cell_tolerance=cell_tolerance,
            kill_timeout=kill_timeout,
            parse_vaspout=parse_vaspout,
            # host=host,
            # unixsocket=unixsocket,
            # port=port,
            # timeout=timeout,
            # log=log,
        )
        input_params.update(**kwargs)

        # Store the original self.command to self._vasp_command since
        # FileIOClientLauncher makes use of it. Leave the port & unixsocket as formatters
        # When starting the socket-I/O interface like SocketIOCalculator(calc=calc, **params)
        # there is no need to specify any socket parameters.
        self._vasp_command = self.command
        self.command = f"{sys.executable} -m vasp_interactive.socketio -p {{port}} -sn {{unixsocket}} -ht {host}"
        self._ensure_directory()
        # save vpi settings
        param_file = self._indir(".vpi_params.pkl")
        with open(param_file, "wb") as f:
            pickle.dump(input_params, f)

        if use_socket:
            self.socket_client = VPISocketClient(
                host=host,
                port=port,
                unixsocket=unixsocket,
                timeout=timeout,
                log=log,
                comm=None,
            )
            self.socket_client.attach_parent_calc(self)
        else:
            self.socket_client = None
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

    def _txt_to_handler(self):
        """Wrap the self.txt into a file handler for parsing
        only relevant for read_results when OUTCAR and vasprun.xml are both trimmed (vasp 5.x)
        """

        txt = self.txt
        if txt is None:
            fd = None
        else:
            if isinstance(txt, str):
                if txt == "-":
                    fd = None
                else:
                    self._ensure_directory()
                    txt = self._indir(txt)
                    fd = open(txt, "r")
            elif hasattr(txt, "read"):
                fd = txt
            else:
                raise RuntimeError(
                    "txt should either be a string"
                    "or an I/O stream, got {}".format(txt)
                )
        return fd

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
            out.flush()

    def _start_vasp_process(self, atoms, out):
        """Helper function to be used inside _run() to make sure the VASP process is running.
        A new VASP process will be started if:
        1. No VASP process has been started
        2. The old VASP process successfully returned but more positions are provided
        If the VASP process is running, do nothing.
        Return of the function ensures self.process is not None
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
            # Dynamic generation of command args
            command = self.make_command(self._vasp_command)
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
            self.pid = self.process.pid
            self.steps = 0
        else:
            retcode = self.process.poll()
            # Still running, do nothing
            if retcode is None:
                return
            if (retcode == 0) and (self.allow_restart_process):
                pid = self.process.pid
                self.process = None
                self.pid = None
                self._stdout(
                    "It seems your VASP process exited normally. I'll restart a new one.",
                    out=out,
                )
                warn(
                    (
                        f"VASP process (pid={pid}) exits normally but new positions are still provided. "
                        "A new VASP process will be started. "
                        "To supress this warning, you may want to increase the NSW number in your settings."
                    )
                )
                # Call the start process again
                self._start_vasp_process(atoms=atoms, out=out)
            else:
                raise RuntimeError(
                    (
                        f"The VASP process has exited with code {retcode} but "
                        "you're still providing new atom positions. "
                        "If the return code is 0, you can try to set allow_restart_process=True "
                        "to enable auto restart of the VASP process."
                    )
                )
        return

    def _write_atoms_stdin(self, atoms, out, require_cell_stdin):
        """Helper function to write input positions / cell in _run().
        This function adapts to different VASP compilations.
        """
        self._stdout("Inputting positions...\n", out=out)
        # ase atoms --> self.sort --> write to position for VASP
        # use scaled positions wrap back to cell
        for atom in atoms.get_scaled_positions()[self.sort]:
            self._stdin(" ".join(map("{:19.16f}".format, atom)), out=None)
        # INPOS subroutine does not split the positions, so manually write into vasp.out
        self._stdout("New scaled positions\n", out=out)
        for i in range(len(atoms)):
            self._stdout(self.process.stdout.readline(), out=out)
        self._stdout("Old scaled positions\n", out=out)
        for i in range(len(atoms)):
            self._stdout(self.process.stdout.readline(), out=out)

        # An additional line "POSITIONS: read from stdin"
        text = self.process.stdout.readline()
        self._stdout(text, out=out)
        assert "POSITIONS: read from stdin" in text

        # Determine if lattice input is needed
        text = self.process.stdout.readline()
        self._stdout(text, out=out)
        if "LATTICE: reading from stdin" in text:
            for vec in atoms.cell:
                self._stdin(" ".join(map("{:19.16f}".format, vec)), out=None)
            # Finish the lattice vector outputs. Note there can be multiple empty lines
            # In total, there should be 9 lines before the closing sentence
            count = 0
            while count < (2 * 3 + 3):
                text = self.process.stdout.readline()
                self._stdout(text, out=out)
                if len(text.strip()) != 0:
                    count += 1
            assert "LATTICE: read from stdin" in text
        elif require_cell_stdin:
            # Cannot continue if cell change is required but VASP does not support
            raise RuntimeError(
                (
                    "The unit cell changes in your calculation but VASP does not support writing lattice parameters to stdin. "
                    "Please consider applying this patch https://github.com/ulissigroup/vasp-interactive first."
                )
            )
        return

    def _run(self, atoms, out, require_cell_stdin):
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
        # 1. Start / renew the VASP process
        self._start_vasp_process(atoms=atoms, out=out)

        # 2. Input the new positions (position and / or cell) if steps > 0
        if self.steps > 0:
            self._write_atoms_stdin(
                atoms=atoms, out=out, require_cell_stdin=require_cell_stdin
            )

        # 3. Start read cycle until the current ionic step
        #    finishes by "POSITIONS: reading from stdin"
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text, out=out)
            # Read vasp version from stdio, warn user of VASP5 issue
            if self._read_vasp_version_stream(text):
                if _int_version(self.version) < 6:
                    warn(
                        (
                            "VASP 5.x. may be not fully compatible with VaspInteractive."
                            "See this issue for details https://github.com/ulissigroup/vasp-interactive/issues/6. "
                            "If you encounter similar error messages, try using VASP 6.x or apply our patch if possible."
                        )
                    )
            if "POSITIONS: reading from stdin" in text:
                return

        # 4. VASP process exits. Check if everything's running ok
        retcode = self.process.poll()
        if retcode == 0:
            self._stdout("VASP terminated normally\n", out=out)
            if self.steps > self.incar_nsw:
                self._stdout(
                    (
                        "However the maximum ionic iterations have been reached. "
                        "Consider increasing NSW number in your calculation.\n"
                    ),
                    out=out,
                )
        else:
            # If we've reached this point, then VASP has exited without asking for
            # new positions, meaning it either exited without error unexpectedly,
            # or it exited with an error. Either way, we need to raise an error.
            raise RuntimeError(
                "VASP exited unexpectedly with exit code {}" "".format(retcode)
            )
        return

    def close(self):
        """Soft stop approach for the stream-based VASP process
        Works by writing STOPCAR file and runs two dummy scf cycles
        Note: if MLFF is enabled, STOPCAR may not work when energy/forces are
        infered from ML model. In that case, simply send SIGINT and wait till timeout
        """
        # Make sure the MPI process is awake before the termination
        if self.pause_mpi:
            # If the mpi_state is never evaluated (None) or RUNNING,
            # we don't need to continue the mpi process
            if getattr(self, "mpi_state", None) == "PAUSED":
                self._resume_calc()

        if self.process is None:
            self.pid = None
            if hasattr(self, "mpi_match"):
                self.mpi_match = None
                self.mpi_state = None
        elif self.process.poll() is not None:
            # For whatever reason the vasp process stops prematurely (possibly too small nsw)
            # do a clean up
            retcode = self.process.poll()
            with self._txt_outstream() as out:
                self._stdout(f"VASP exited with code {retcode}.", out=out)
            self._reset_process()
        elif self._use_mlff():
            # Send only SIGINT to process to terminate.
            # Current with VASP 6.3.0
            warn("Terminate MLFF-VASP with SIGINT. OUTCAR may be incomplete.")
            self._send_mpi_signal(signal.SIGINT)
            self._reset_process()
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
                    self._run(self.atoms, out=out, require_cell_stdin=False)

                    if self.process.poll() is not None:
                        self._stdout(
                            f"VASP exited with code {self.process.poll()}.", out=out
                        )
                        self._reset_process()
                        return
                while self.process.poll() is None:
                    time.sleep(1)
                self._stdout("VASP has been closed\n", out=out)
                self._reset_process()
        return

    def _send_mpi_signal(self, sig):
        """Send signal to the mpi process within self.process
        If the process cannot be found, return without affecting the state
        """
        # Whenever cannot locate the pid via psutil, reset the VASP process
        try:
            pid = self.process.pid
            psutil_proc = psutil.Process(pid)
        except Exception as e:
            warn("VASP process no longer exists. Will reset the calculator.")
            self._reset_process()
            return

        if (self.pid == pid) and getattr(self, "mpi_match", None) is not None:
            match = self.mpi_match
        else:
            self.pid = pid
            match = _find_mpi_process(pid)
            self.mpi_match = match
        if (match["type"] is None) or (match["process"] is None):
            warn(
                "Cannot find the mpi process or you're using different ompi wrapper. Will not send stop signal to mpi."
            )
            return
        elif match["type"] == "mpi":
            mpi_process = match["process"]
            mpi_process.send_signal(sig)
        elif match["type"] == "slurm":
            slurm_step = match["process"]
            _slurm_signal(slurm_step, sig)
        else:
            raise ValueError("Unsupported process type!")
        return

    def _pause_calc(self, sig=signal.SIGTSTP):
        """Pause the vasp processes by sending SIGTSTP to the master mpirun process
        If the current pid are the same with previous record, do not query the mpi pid or slurm stepid
        """
        # Always check if current calculator allows pause in case external optimizers
        # use self._resume_calc / self._pause_calc pair without the context
        if not self.pause_mpi:
            return

        if not self.process:
            return

        self._send_mpi_signal(sig=sig)
        self.mpi_state = "PAUSED"
        return

    def _resume_calc(self, sig=signal.SIGCONT):
        """Resume the vasp processes by sending SIGCONT to the master mpirun process
        If the current pid are the same with previous record, do not query the mpi pid or slurm stepid
        """
        # Always check if current calculator allows pause in case external optimizers
        # use self._resume_calc / self._pause_calc pair without the context
        if not self.pause_mpi:
            return

        if not self.process:
            return

        self._send_mpi_signal(sig=sig)
        self.mpi_state = "RUNNING"
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
        system_changes=all_changes,
    ):
        check_atoms(atoms)

        # TODO: check if this is an old hack left from ALMLP days?
        if hasattr(self, "system_changes") and self.system_changes is not None:
            system_changes = self.system_changes
            self.system_changes = None

        if not system_changes:
            # No need to calculate, calculator has stored calculated results
            return

        # VaspInteractive supports positions change by default.
        # Change of cell is supported by a patch provided by this package, but needs to check on the fly
        # TODO: send flag to _run to warn cell change
        if "numbers" in system_changes:
            if self.process is not None:
                raise NotImplementedError(
                    (
                        "VaspInteractive does not support change of chemical formula. "
                        "Please create a new calculator instance or use standard Vasp calculator"
                    )
                )

        # Does the VASP interactive mode needs to take cell as inputs?
        # If require_cell_stdin is True, _run must get a "LATTICE: reading from stdin" line
        # after sending in the positions, otherwise throws NotImplementedError
        # Otherwise, sending cell to stdin is optional
        if ("cell" in system_changes) and (self.process is not None):
            require_cell_stdin = True
        else:
            require_cell_stdin = False
        # elif "cell" in system_changes:

        # if self.process is not None:
        #     raise NotImplementedError(
        #         (
        #             "VaspInteractive does not support change of lattice parameters. "
        #             "Set VaspInteractive.cell_tolerance to a higher value if you think it's caused by round-off error. "
        #             "Otherwise, please create a new calculator instance or use standard Vasp calculator"
        #         )
        #     )

        self.clear_results()
        if atoms is not None:
            self.atoms = atoms.copy()

        with self._txt_outstream() as out:
            self._run(self.atoms, out=out, require_cell_stdin=require_cell_stdin)
            self.steps += 1
            # special condition: job runs with nsw limit reached.
            # In the interactive mode, VASP won't exit until steps > nsw
            # this can be problematic if the next user input is a different position
            # so we simply run a dummy step using current positions to gracefully terminate VASP
            if self.steps >= self.incar_nsw:
                self._run(self.atoms, out=out, require_cell_stdin=require_cell_stdin)

        # Use overwritten `read_results` method
        self.read_results()
        return

    def read_results(self):
        """Overwrites the `read_results` from parent class.
        Depending on the state of the VASP patch, the vasprun.xml, OUTCAR or vasp.out fields
        may or may not be truncated after each ionic step. `read_results` will first try
        xml parsing using parent class's read_results.
        If not working then OUTCAR and finally read from vasp.out if all above failed.

        The xml calc reader self._xml_calc may be generated by a place-holder SP calculator
        to make _set_old_parameters working.
        """
        # Special notes for VASP 6.3+ with MLFF flag on. Currently ase does not
        # handle VASP output with MLFF settings, so we always need to parse from OUTCAR
        # with some pre-filtration of text.
        # Good news is, the OUTCAR files are always completed :)
        if self._use_mlff():
            use_mlff = True
            self._xml_complete = False
            warn(
                "Parsing output from MLFF-enabled VASP is experimental. Use at your own risk!"
            )
        else:
            use_mlff = False
        # Record if xml contents are complete at each ionic step
        if self._xml_complete is not False:
            try:
                super().read_results()
                self._xml_complete = True
            except Exception as e:
                warn(
                    (
                        "Direct reading from xml files failed. "
                        "You VASP build may be truncating outputs."
                        "I'll recover the remaining information using OUTCAR and vasp.out at this point"
                    )
                )
                self._xml_complete = False

        # Some properties, like E-Fermi might not exist in vasprun.xml but rather in OUTCAR
        # Try to complete as much as possible
        outcar = self.load_file("OUTCAR")
        if use_mlff:
            # breakpoint()
            outcar = _preprocess_mlff_outcar(outcar)
        properties = ["stress", "fermi", "nbands", "dipole"]
        for prop in properties:
            if prop not in self.results.keys():
                try:
                    # use read_xxx method to parse outcar
                    result = getattr(self, f"read_{prop}")(lines=outcar)
                    self.results[prop] = result
                except Exception:
                    pass
        # TODO: fermi should be added to _xml_calc instead
        # breakpoint()
        # No need to further parsing if xml file parsing finished
        if self._xml_complete:
            return

        self.converged = self.read_convergence(lines=outcar)
        self.version = self.read_version()

        # Do a quick test to see if OUTCAR is truncated
        # remember the status of _outcar_complete later
        if self._outcar_complete is None:
            e, fe = self.read_energy(lines=outcar)
            f = self.read_forces(lines=outcar)
            if (e == 0) or (fe == 0) or (f is None):
                self._outcar_complete = False
                warn(
                    (
                        "OUTCAR contents are truncated and I'll try to parse "
                        "energy and forces from vasp.out. Note some information "
                        "like stress may be incorrect."
                    )
                )
            else:
                self._outcar_complete = True

        if (self._outcar_complete is False) and (
            getattr(self, "parse_vaspout", False) is False
        ):
            raise RuntimeError(
                (
                    "Cannot continue because both vasprun.xml and OUTCAR are truncated. "
                    "If you are running with VASP 5.x, set parse_vaspout=True. "
                )
            )

        # Energy and magmom have multiple return values
        if "free_energy" not in self.results.keys():
            try:
                energy_free, energy_zero = self.read_energy(lines=outcar, all=False)
                self.results.update(dict(free_energy=energy_free, energy=energy_zero))
            except Exception as e:
                raise RuntimeError(("Failed to obtain energy from calculator.")) from e

        if "forces" not in self.results.keys():
            try:
                # self.read_forces already contain sort information
                forces = self.read_forces(lines=outcar)
                self.results.update(dict(forces=forces))
            except Exception as e:
                raise RuntimeError(("Failed to obtain forces from calculator.")) from e

        if "magmom" not in self.results.keys():
            try:
                magmom, magmoms = self.read_mag(lines=outcar)
                self.results.update(dict(magmom=magmom, magmoms=magmoms))
            except Exception:
                pass

        # print(self.results)
        # print(all_properties)
        # Construct a fake _xml_calc object
        # results = self.results.copy()
        # nbands = results.pop("nbands", None)
        # print(self.results["forces"])
        results = {k: v for k, v in self.results.items() if k in all_properties}
        # print(list(results.keys()))
        calc_xml = SinglePointDFTCalculator(
            atoms=self.atoms, efermi=self.results.get("fermi", None), **results
        )
        self._xml_calc = calc_xml
        # Manunal old keywords handling
        self._set_old_keywords()

        # Store the parameters used for this calculation
        self._store_param_state()

    def read_energy(self, all=False, lines=None):
        """Overwrite the parent read_energy
        VASP 5.x output behavior is unpredictable and should always use the value in vasp.out
        parameter vasp5 enforces read using vasp.out (or any txt) and uses only all=False
        """
        if self._outcar_complete is not False:
            try:
                fe, e0 = super().read_energy(all=all, lines=lines)
            except Exception:
                fe, e0 = [0, 0]
        else:
            if getattr(self, "parse_vaspout", False) is False:
                raise RuntimeError(
                    (
                        "Cannot parse potential energy. Most likely vasprun.xml and OUTCAR are both corrupt. \n"
                        "If you are running with VASP 5.x, add parse_vaspout=True to the calculator "
                    )
                )
            f_vaspout = self._txt_to_handler()
            if f_vaspout is None:
                raise RuntimeError(
                    (
                        "Fails to locate file for vasp output. Please don't set the calculator's txt parameter to either '-' or None"
                    )
                )
            vaspout = f_vaspout.readlines()
            f_vaspout.close()
            try:
                fe, e0 = parse_vaspout_energy(vaspout, all=False)
            except Exception as e:
                raise RuntimeError(("Cannot parse energy from vasp output.")) from e
        return fe, e0

    def read_forces(self, all=False, lines=None, vasp5=False):
        """Overwrite the parent read_forces
        VASP 5.x output behavior is unpredictable and should always use the value in vasp.out
        parameter vasp5 enforces read using vasp.out (or any txt)
        """
        if self._outcar_complete is not False:
            try:
                forces = super().read_forces(all=all, lines=lines)
            except Exception:
                forces = None

        # Upstream read_energy from ase has a flaw that returns 0, 0 if parsing failed
        # Need to handle such case
        else:
            if getattr(self, "parse_vaspout", False) is False:
                raise RuntimeError(
                    (
                        "Cannot parse forces. Most likely vasprun.xml and OUTCAR are both corrupt. \n"
                        "If you are running with VASP 5.x, add parse_vaspout=True to the calculator "
                    )
                )
            f_vaspout = self._txt_to_handler()
            if f_vaspout is None:
                raise RuntimeError(
                    (
                        "Fails to locate file for vasp output. Please don't set the calculator's txt parameter to either '-' or None"
                    )
                )
            vaspout = f_vaspout.readlines()
            f_vaspout.close()
            try:
                # Remember to resort the forces according to ase input
                forces = parse_vaspout_forces(vaspout, all=False)[self.resort]
            except Exception as e:
                raise RuntimeError(("Cannot parse forces from vasp output.")) from e
        return forces

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
        if (self.socket_client is not None) and (not self.socket_client.closed):
            self.socket_client.close()
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
                    f"{traceback.format_exc(limit=2)}\n"
                    "Will now force closing the VASP process. "
                    "The OUTCAR and vasprun.xml outputs may be incomplete"
                ),
                file=sys.stderr,
            )
            if self.process is not None:
                if self.process.poll() is None:
                    # self.process.kill()
                    self._send_mpi_signal(signal.SIGKILL)
            # Reset process tracker
            self._reset_process()
        return

    def _reset_process(self):
        """Reset the record for process in the calculator.
        Useful if the process is missing or reset the calculator.
        """
        # Reset process tracker
        self.process = None
        self.pid = None
        if hasattr(self, "mpi_match"):
            self.mpi_match = None
            self.mpi_state = None

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
        new._reset_process()
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

    def _use_mlff(self):
        """Parameter check if mlff is used in this calculation.
        ml_lmlff is currently not recognized as a valid parameter so need to search
        the "custom" key.
        """
        params = self.parameters.copy()
        # Reserved for future releases of ase
        if params.get("ml_lmlff", False) is True:
            return True

        if "custom" not in params:
            return False

        return params["custom"].get("ml_lmlff", False)

    def are_results_dft(self, lines=None):
        """Parse the OUTCAR file to see if last returned energy/forces are from DFT or ML inference
        This function can be called at the end of each ionic step
        """
        if not self._use_mlff():
            return True
        if not lines:
            lines = self.load_file("OUTCAR")
        pat_all = r"energy\s+without\s+entropy"
        pat_ml = r"^\s+ML\s+energy\s+without\s+entropy"
        last_line = None
        for line in lines:
            if re.search(pat_all, line):
                last_line = line
        if last_line and re.search(pat_ml, last_line):
            return False
        else:
            return True


    # socket-related
    def irun(self, atoms, use_stress=None):
        """Make the client run in iterative mode"""
        if self.socket_client is None:
            raise NotImplementedError(
                "Cannot use socket io mode without specifying use_socket=True"
            )
        atoms.calc = self
        return self.socket_client.irun(atoms, use_stress=use_stress)

    def run(self, atoms, use_stress=None):
        """Infinitely run client code"""
        if self.socket_client is None:
            raise NotImplementedError(
                "Cannot use socket io mode without specifying use_socket=True"
            )
        atoms.calc = self
        # breakpoint()
        self.socket_client.run(atoms, use_stress=use_stress)
        return
