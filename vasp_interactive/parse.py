"""Functions parsing OUTCAR or vasp.out files which are not present in parent
Vasp calculator but can he helpful for job diagnosis
"""
import time
import os
import sys
import psutil
import signal
import re
import numpy as np


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


def parse_vaspout_energy(lines, all=False):
    """Parse the energy and force lines from lines or stdout
    This should only be relevant when vasp5.x is involved, as a last resort
    when both OUTCAR and vaspun.xml are truncated. Please use with care.

    The vaspout format are almost identical to OSZICAR https://www.vasp.at/wiki/index.php/OSZICAR
    FORCES:
        (N x 3) fields of force
    N F= XX E0= XX  d E =XX
    The F & E0 line should most likely maintain the same format but there can be extra output lines like
    vdW dispersion etc. Some codes can be traken from https://github.com/materialsproject/pymatgen/blob/v2022.9.8/pymatgen/io/vasp/outputs.py#L4253-L4394

    For upstream methods in the calculator class, one must first check if the calc.txt is neither "-" nor None,
    otherwise the vaspout lines cannot be parsed.
    """
    free_energies, zero_energies = [], []
    ionic_pattern = re.compile(
        (
            r"(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
            r"E0=\s*([\d\-\.E\+]+)\s+"
            r"d\s*E\s*=\s*([\d\-\.E\+]+)$"
        )
    )
    for line in lines:
        line = line.strip()
        if ionic_pattern.match(line.strip()):
            m = ionic_pattern.match(line.strip())
            free_energies.append(float(m.group(2)))
            zero_energies.append(float(m.group(3)))
    if all is False:
        free_energies = free_energies[-1]
        zero_energies = zero_energies[-1]
    return free_energies, zero_energies


def parse_vaspout_forces(lines, all=False):
    """Parse the energy and force lines from lines or stdout
    This should only be relevant when vasp5.x is involved, as a last resort
    when both OUTCAR and vaspun.xml are truncated. Please use with care.

    The vaspout format are almost identical to OSZICAR https://www.vasp.at/wiki/index.php/OSZICAR
    FORCES:
        (N x 3) fields of force
    N F= XX E0= XX  d E =XX
    The F & E0 line should most likely maintain the same format but there can be extra output lines like
    vdW dispersion etc. Some codes can be taken from https://github.com/materialsproject/pymatgen/blob/v2022.9.8/pymatgen/io/vasp/outputs.py#L4253-L4394

    For upstream methods in the calculator class, one must first check if the calc.txt is neither "-" nor None,
    otherwise the vaspout lines cannot be parsed.
    """
    forces = []
    current_line = 0
    for i, line in enumerate(lines):
        if i < current_line:
            continue
        else:
            current_line = i
        if "FORCES:" in line:
            forces_this_step = []
            while True:
                current_line += 1
                try:
                    fx, fy, fz = np.fromstring(
                        lines[current_line], dtype=float, sep=" "
                    )
                    forces_this_step.append([fx, fy, fz])
                except Exception:
                    break
            forces.append(forces_this_step)
    # forces = np.array(forces)
    if all is False:
        forces = np.array(forces[-1])
    else:
        # May be unusable if vasp.out is appended from previous calculation.
        # Use all=True with caution
        forces = np.array(forces)
    return forces
