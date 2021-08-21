"""Testing if environmental variable correctly passes to process
"""
import pytest
from vasp_interactive import VaspInteractive
import tempfile
from pathlib import Path
import os
import sys


def test_command_env():
    from ase.build import molecule

    # When testing under the github action system, it should be
    # `mpirun -np 1 <options> <path>/vasp_gam`
    default_command = os.environ["VASP_COMMAND"]

    """Test if VaspInteractive correctly catches command if env changes
    """
    atoms = molecule("H2", vacuum=5, pbc=True)
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        calc = VaspInteractive(xc="pbe", directory=tempdir)
        with calc:
            atoms.calc = calc
            assert calc.command is None
            command = calc.make_command()
            assert command == default_command
            atoms.get_potential_energy()
            # calc.process the shell mode, so args is a string
            assert calc.process.args == default_command

            # Harmless change
            new_command = default_command.replace("-np", "-n")
            calc.command = new_command
            # Should be no change
            assert calc.process.args == default_command

        with calc:
            atoms.rattle()
            atoms.get_potential_energy()
            print(calc.process.args)
            assert calc.process.args == new_command

        with calc:
            # reset the environ
            os.environ["VASP_COMMAND"] = "echo TEST && " + default_command
            atoms.rattle()
            calc.command = None
            atoms.get_potential_energy()
            assert calc.command is None
            print(calc.process.args)
            assert "TEST" in calc.process.args
