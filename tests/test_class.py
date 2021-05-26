import pytest
from vasp_interactive import VaspInteractive
import tempfile
from pathlib import Path

def test_class():
    from ase.build import molecule
    """Test if VaspInteractive correctly write inputs
    """
    atoms = molecule("C2H2", vacuum=5)
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        calc = VaspInteractive(xc="pbe", 
                               directory=tempdir)
        # Initialization
        assert calc.process is None
        assert hasattr(calc.txt, "write")
        assert calc.txt.mode == "a"
        assert calc.txt.closed is False
        
        # Write inputs
        calc.write_input(atoms)
        incar = tempdir / "INCAR"
        assert incar.is_file()
        with open(incar, "r") as fd:
            context = fd.read()
            print(context)
        assert "INTERACTIVE = .TRUE." in context
    return
    