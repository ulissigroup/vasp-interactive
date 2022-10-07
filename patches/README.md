# Patch the VASP source code for better interactive behavior

------

**!!!DANGER ZONE!!!**

The tutorial in this readme file will change the official VASP source code. Although 
these patches are meant only to modify the behavior of the interactive mode, sudden change / incorrect 
formatting may cause build to fail or unexpected program behavior. Please check the contents of the patches 
before applying them.

------

## Do I need these patches?

If you're only using `VaspInteractive` for relaxation tasks which does not involve lattice change 
(e.g. slab geometry optimization, molecular dynamics, etc.), compiling VASP from its default source code is sufficient. 
However, you may consider build a patched VASP version for the following reasons:

1. No unexpected truncation of output file contents (especially for VASP 5.x)
2. Supporting input of lattice change (e.g. performing lattice optimization, fitting EOS curve)
3. Full compatibility with the iPI socket interface via `VaspInteractive`

Note you can also compile the VASP source code to support direct iPI protocol using patches provided with the 
[iPI package](https://github.com/i-pi/i-pi/tree/master/examples/VASP) but currently limited to VASP 5.x.

## Howto

First, locate the directory of uncompressed original VASP source code, for example `vasp.6.3.0/src`, 
which should contain the `.F` Fortran files.
Simply use `patch.py` to apply the changes
```bash
git clone https://github.com/ulissigroup/vasp-interactive.git
cd vasp-interactive/patches
# Change the directory accordingly
python patch.py vasp.6.3.0/src
```

Confirm the new subroutine `INLATT` is patched correctly:
```bash
grep INLATT -5 src/*.F
```

Then compile VASP using the normal approach, such as:
```bash
make -j4 all
```

## Version compatibility

The section for interactive mode has not been significantly changed in VASP source code since 
version 5.4, meaning the patches may be valid for several future releases. Currently tested versions are:
- 5.4.4pl2
- 6.1.2
- 6.2.0
- 6.3.0

Please contact the maintainer if you have issues applying the patches due to VASP source code change.




