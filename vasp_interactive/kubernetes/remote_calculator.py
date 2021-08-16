"""Running vasp interactive in separate remote resources, other than the local processes
   Currently the following remote(s) are implemented
   - kubernetes
   
   Be careful when using this module for production purpose.
"""
from vasp_interactive import VaspInteractive

class KubeVaspInteractive(VaspInteractive):
    pass