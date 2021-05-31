try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name="vasp-interactive",
    version="0.0.1",
    packages=[
        "vasp_interactive",
    ],
    install_requires=[
        "ase",
        # "pymatgen",
    ],
)
