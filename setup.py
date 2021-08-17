try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name="vasp-interactive",
    version="0.0.1",
    packages=["vasp_interactive", "vasp_interactive.kubernetes"],
    install_requires=[
        "ase",
        # "pymatgen",
    ],
)
