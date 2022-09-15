import time
import os
import sys
import psutil
import signal
import re
import numpy as np
from contextlib import contextmanager

# Timeout in seconds before a "forced" kill
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

    If srun is found as the process, need to use `scancel` to pause / resume the job step
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
            match["process"] = proc
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
        match["process"] = mpi_proc = mpi_candidates[-1]

    return match


def _int_version(version_string):
    """Convert version string to its major version
    5.4.4pl2 --> 5
    6.3.0-gpu --> 6
    """
    major = int(version_string.split(".")[0])
    return major
