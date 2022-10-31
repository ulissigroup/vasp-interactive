"""Testing if VaspInteractive correctly handles symmetry change
"""
import pytest
from ase.build import molecule
from ase.calculators.vasp import Vasp
from vasp_interactive import VaspInteractive
import tempfile
import os
import sys
from _common_utils import skip_lattice_if_incompatible

h2_origin = molecule("H2", cell=[8, 8, 8], pbc=True)
# Slightly change cell so system think h2_1 is still the same as h2_origin
# only in VaspInteractive. Will become different in default Vasp
h2_1 = molecule("H2", cell=[8, 8, 8 + 1e-10], pbc=True)
# Change cell but h2_2 is different
h2_2 = molecule("H2", cell=[8, 8, 8 + 1e-7], pbc=True)

# all default settings
vasp_params = dict(xc="pbe", kpts=(1, 1, 1), gamma=True)


def test_check_state():
    """Unit test for check_state function. Only test when no lattice info needed."""
    skip_lattice_if_incompatible(reverse=True)
    vasp = Vasp(**vasp_params)
    vpi = VaspInteractive(**vasp_params)

    vasp.atoms = h2_origin
    # Both will report a cell change
    system_changes = vasp.check_state(h2_1)
    assert "positions" not in system_changes
    assert "cell" in system_changes

    system_changes = vasp.check_state(h2_2)
    assert "positions" not in system_changes
    assert "cell" in system_changes

    # for VaspInteractive, h2_1 is accepted
    vpi.atoms = h2_origin
    system_changes = vpi.check_state(h2_1)
    assert "positions" not in system_changes
    assert "cell" not in system_changes

    system_changes = vpi.check_state(h2_2)
    assert "positions" not in system_changes
    assert "cell" in system_changes


def test_calculation():
    skip_lattice_if_incompatible(reverse=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        vpi = VaspInteractive(directory=tmpdir, **vasp_params)
        with vpi:
            h2_origin.calc = vpi
            e1 = h2_origin.get_potential_energy()
            h2_1.calc = vpi
            try:
                e2 = h2_1.get_potential_energy()
            except Exception:
                pytest.fail("h2_1 cell tolerance should be accepted")
            h2_2.calc = vpi
            with pytest.raises(Exception):
                e3 = h2_2.get_potential_energy()

            assert pytest.approx(e1) == e2


def test_always_true():
    pass
