import pytest
from vasp_interactive import VaspInteractive
import tempfile
from pathlib import Path
import os
import sys


def test_class():
    from ase.build import molecule

    """Test if VaspInteractive correctly write inputs
    """
    atoms = molecule("C2H2", vacuum=5, pbc=True)
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        calc = VaspInteractive(xc="pbe", directory=tempdir)
        # Initialization
        assert calc.process is None

        # txt io
        with calc._txt_outstream() as out:
            assert out.mode == "a"
            assert hasattr(out, "write")
            assert os.path.exists(out.name)
            print(out.name)

        # Write inputs
        calc.write_input(atoms)
        incar = tempdir / "INCAR"
        assert incar.is_file()
        with open(incar, "r") as fd:
            context = fd.read()
            print(context)

        # Additional parameters for VaspInteractive
        assert "INTERACTIVE = .TRUE." in context
        assert "IBRION = -1" in context
        assert "POTIM = 0" in context
        assert "IWAVPR = 11" in context

        # Try manual closing, since no process is given
        # Should we close txt at the same time?
        assert calc.close() is None
        #         assert calc.txt.closed is False

        # Context manager will close txt correctly
        with calc:
            pass
        assert calc.process is None
    return


def test_output():
    from ase.build import molecule

    """Test if VaspInteractive correctly write inputs
    """
    atoms = molecule("C2H2", vacuum=5, pbc=True)
    with tempfile.TemporaryDirectory() as tempdir:
        # No value provided, default to vasp.out
        calc = VaspInteractive(xc="pbe", directory=tempdir)
        with calc._txt_outstream() as out:
            assert "vasp.out" in out.name

        # use std
        calc = VaspInteractive(xc="pbe", directory=tempdir, txt="-")
        with calc._txt_outstream() as out:
            assert out == sys.stdout

        # use no output
        calc = VaspInteractive(xc="pbe", directory=tempdir, txt=None)
        with calc._txt_outstream() as out:
            assert out is None
