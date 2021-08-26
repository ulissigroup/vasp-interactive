"""Testing if VaspInteractive correctly handles symmetry change
"""
import pytest
from ase.build import graphene
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
import tempfile
import os
import sys

gr = graphene(vacuum=10)
gr.pbc = True
# drift the carbon atom just to perturb the symmetry
gr_distorted = gr.copy()
gr_distorted[-1].x += 0.5

# all default settings
vasp_params = dict(xc="pbe", kpts=(5, 5, 1), gamma=True)

# stop criteria
stop_crit = 5e-3

with tempfile.TemporaryDirectory() as tmpdir:
    calc_vasp = Vasp(directory=tmpdir, **vasp_params)
    e_gr0 = calc_vasp.get_potential_energy(gr.copy())
    e_gr1 = calc_vasp.get_potential_energy(gr_distorted.copy())


def test_sym_on():
    """Test if symmetry is set to on"""
    with tempfile.TemporaryDirectory() as tmpdir:
        calc_vpi = VaspInteractive(directory="./", isym=2, **vasp_params)
        # First run wavefunction on high sym and then low sym
        with calc_vpi:
            e0 = calc_vpi.get_potential_energy(gr.copy())
            e1 = calc_vpi.get_potential_energy(gr_distorted.copy())
            print(e_gr0, e0)
            print(e_gr1, e1)
            assert pytest.approx(e_gr0, stop_crit) == e0
            assert pytest.approx(e_gr1, stop_crit) != e1

        # First run wavefunction on low sym and then high sym
        with calc_vpi:
            e1 = calc_vpi.get_potential_energy(gr_distorted.copy())
            e0 = calc_vpi.get_potential_energy(gr.copy())
            print(e_gr0, e0)
            print(e_gr1, e1)
            assert pytest.approx(e_gr0, stop_crit) != e0
            assert pytest.approx(e_gr1, stop_crit) == e1


def test_sym_off():
    """Test if symmetry is set to on"""
    with tempfile.TemporaryDirectory() as tmpdir:
        calc_vpi = VaspInteractive(directory="./", **vasp_params)
        # First run wavefunction on high sym and then low sym
        with calc_vpi:
            e0 = calc_vpi.get_potential_energy(gr.copy())
            e1 = calc_vpi.get_potential_energy(gr_distorted.copy())
            print(e_gr0, e0)
            print(e_gr1, e1)
            assert pytest.approx(e_gr0, stop_crit) == e0
            assert pytest.approx(e_gr1, stop_crit) == e1

        # First run wavefunction on low sym and then high sym
        with calc_vpi:
            e1 = calc_vpi.get_potential_energy(gr_distorted.copy())
            e0 = calc_vpi.get_potential_energy(gr.copy())
            print(e_gr0, e0)
            print(e_gr1, e1)
            assert pytest.approx(e_gr0, stop_crit) == e0
            assert pytest.approx(e_gr1, stop_crit) == e1
