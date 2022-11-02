try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name="vasp-interactive",
    version="0.1.0rc0",
    packages=["vasp_interactive", "vasp_interactive.kubernetes"],
    install_requires=[
        "ase",
        "psutil",
        # "pymatgen",
    ],
)
