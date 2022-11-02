import time
import os
import sys
import psutil
import signal
import re
import numpy as np
from warnings import warn
import subprocess
from contextlib import contextmanager

# Timeout in seconds before a "forced" kill
DEFAULT_KILL_TIMEOUT = 60


def _run_process(commands, shell=False, print_cmd=True, cwd=".", capture_output=False):
    """Wrap around subprocess.run
    Returns the process object
    """
    full_cmd = " ".join(commands)
    if print_cmd:
        print(" ".join(commands))
    if shell is False:
        proc = subprocess.run(
            commands, shell=shell, cwd=cwd, capture_output=capture_output
        )
    else:
        proc = subprocess.run(
            full_cmd, shell=shell, cwd=cwd, capture_output=capture_output
        )
    if proc.returncode == 0:
        return proc
    else:
        raise RuntimeError(f"Running {full_cmd} returned error code {proc.returncode}")


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

    If srun is found as the process, need to use `scancel` to pause / resume the job step
    """
    allowed_names = set(["mpirun", "mpiexec", "orterun", "oshrun", "shmemrun"])
    allowed_vasp_names = set(["vasp_std", "vasp_gam", "vasp_ncl"])
    if mpi_program:
        allowed_names.add(mpi_program)
    if vasp_program:
        allowed_vasp_names.add(vasp_program)
    try:
        process_list = [psutil.Process(pid)]
    except psutil.NoSuchProcess:
        warn("Psutil cannot locate the pid. Your VASP program seems already exited.")
        match = {"type": None, "process": None}
        return match

    process_list.extend(process_list[0].children(recursive=True))
    mpi_candidates = []
    match = {"type": None, "process": None}
    for proc in process_list:
        name = proc.name()
        if name in ["srun"]:
            # TODO: after the slurm parts get mature the warning can be removed
            warn(
                "VASP is started by slurm's srun command which requires different pause / resume mechanism. "
                "If you spot weird behaviors, disable the pausing of calculator."
            )
            match["type"] = "slurm"
            match["process"] = _locate_slurm_step(vasp_program=vasp_program)
            break
        elif proc.name() in allowed_names:
            # is the mpi process's direct children are vasp_std?
            children = proc.children()
            if len(children) > 0:
                if children[0].name() in allowed_vasp_names:
                    mpi_candidates.append(proc)
    if len(mpi_candidates) > 1:
        warn(
            "More than 1 mpi processes are created. This may be a bug. I'll use the last one"
        )
    if len(mpi_candidates) > 0:
        match["type"] = "mpi"
        match["process"] = mpi_candidates[-1]

    return match


def _int_version(version_string):
    """Convert version string to its major version
    5.4.4pl2 --> 5
    6.3.0-gpu --> 6
    """
    major = int(version_string.split(".")[0])
    return major


def _get_slurm_jobid():
    jobid = os.environ.get("SLURM_JOB_ID", None)
    if jobid is None:
        jobid = os.environ.get("SLURM_JOBID", None)
    return jobid


def _locate_slurm_step(vasp_program="vasp_std"):
    """If slurm job system is involved, search for the slurm step id
    that matches vasp_std (or other vasp commands)

    Steps:
    1. Look for SLURM_JOB_ID in current env
    2. Use `squeue` to locate the vasp_std step (latest)

    squeue
    """
    allowed_vasp_names = set(["vasp_std", "vasp_gam", "vasp_ncl"])
    if vasp_program:
        allowed_vasp_names.add(vasp_program)
    jobid = _get_slurm_jobid()
    if jobid is None:
        # TODO: convert warn to logger
        warn(("Cannot locate the SLURM job id."))
        return None
    # Only 2 column output (jobid and jobname)
    cmds = ["squeue", "-s", "--job", str(jobid), "-o", "%.30i %.30j"]
    proc = _run_process(cmds, capture_output=True)
    output = proc.stdout.decode("utf8").split("\n")
    # print(output)
    candidates = []
    # breakpoint()
    for line in output[1:]:
        try:
            stepid, name = line.strip().split()
        except Exception:
            continue
        if any([v in name for v in allowed_vasp_names]):
            candidates.append(stepid)

    if len(candidates) > 1:
        warn("More than 1 slurm steps are found. I'll use the most recent one")
    if len(candidates) > 0:
        proc = candidates[0]
    else:
        proc = None
    return proc


def _slurm_signal(stepid, sig=signal.SIGTSTP):
    if isinstance(sig, (str,)):
        sig = str(sig)
    elif isinstance(sig, (int,)):
        sig = signal.Signals(sig).name
    else:
        sig = sig.name
    cmds = ["scancel", "-s", sig, str(stepid)]
    proc = _run_process(cmds, capture_output=True)
    output = proc.stdout.decode("utf8").split("\n")
    return

def _preprocess_mlff_outcar(outcar):
    """Pre-filter ML-related lines to normal VASP OUTCAR lines
    This should only work for Vasp.read_X (X=property) method while we don't tackle the 
    read_vasp_out method in ase.io.
    Include (watch for the extra spaces!)
    "^ ML energy  without entropy" --> "^  energy  without entropy"
    "^ free  energy ML TOTEN" --> '^  free  energy   TOTEN'
    "ML energy(sigma->0)" --> "  energy(sigma->0)"
    "TOTAL-FORCE (eV/Angst) (ML)" --> "TOTAL-FORCE (eV/Angst)"

    outcar: full content of outcar
    """
    replacement = [
        (r"^\s+ML\s+energy\s+without\s+entropy", "  energy  without entropy"),
        (r"^\s+free\s+energy\s+ML\s+TOTEN", "  free  energy   TOTEN"),
        (r"ML\s+energy\(sigma->0\)", "   energy(sigma->0)"),
        (r"TOTAL-FORCE\s+\(eV/Angst\)\s+\(ML\)", "TOTAL-FORCE (eV/Angst)"),
        (r"^\s+ML\s+FREE\s+ENERGIE\s+OF\s+THE\s+ION-ELECTRON\s+SYSTEM\s+\(eV\)",
        "  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)"),
        (r"ML\s+FORCE\s+on\s+cell",
        "   FORCE on cell")
    ]
    new_outcar = []
    for line in outcar:
        newline = line
        for pat, rep in replacement:
            newline = re.sub(pat, rep, newline, 0, re.MULTILINE)
        new_outcar.append(newline)
    return new_outcar
