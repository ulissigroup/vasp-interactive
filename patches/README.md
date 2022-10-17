# Patch the VASP source code for better interactive behavior

------

**!!!DANGER ZONE!!!**

You are going to change the official VASP source code. Although 
these patches are meant only to modify the behavior of the interactive mode, sudden change / incorrect 
formatting may cause build to fail or unexpected program behavior. Please check the contents of the patches 
before applying them.

------

## TL; DR

Use `compile_vasp.sh` for the patching and compilation the VASP source code

```
INTERACTIVE_PATCH=patch.py ./compile_vasp.sh <path-to-VASP-source-tgz> <path-to-makefile> <path-to-makefile.include>
```

We provide examples of  `makefile` and `makefile.include` for compilation using PGI compilers in the Nvidia HPCSDK toolchain. 
You can simply adjust the variables to your requirements. Here are a few examples:

- **Example 1**: building pristine VASP 6.3.0 on an 8-core machine, uncompress to `/tmp`, and move the binary files to `/opt/vasp/bin`

```bash
NCORES=8 ROOT=/tmp VASP_BINARY_PATH=/opt/vasp/bin
./compile_vasp.sh vasp.6.3.0.tgz examples/makefile examples/makefile.include.vasp6
```

- **Example 2**: add the interactive mode patch to Example 1

```bash
NCORES=8 ROOT=/tmp VASP_BINARY_PATH=/opt/vasp/bin \
INTERACTIVE_PATCH=patch.py \
./compile_vasp.sh vasp.6.3.0.tgz examples/makefile examples/makefile.include.vasp6
```

You can use any `makefile` and `makefile.include` for your specific compiler-OS combinations.
For a detailed explanation please see the [**How-to**](#howto) section.


## Do I need these patches?

If you're only using `VaspInteractive` for relaxation tasks which does not involve lattice change 
(e.g. slab geometry optimization, molecular dynamics, etc.), compiling VASP from its default source code is sufficient. 
However, you may consider build a patched VASP version for the following reasons:

1. No unexpected truncation of output file contents (especially for VASP 5.x)
2. Supporting input of lattice change (e.g. performing lattice optimization, fitting EOS curve)
3. Full compatibility with the iPI socket interface via `VaspInteractive`

Note you can also compile the VASP source code to support direct iPI protocol using patches provided with the 
[iPI package](https://github.com/i-pi/i-pi/tree/master/examples/VASP) but currently limited to VASP 5.x. Our patch 
focuses only to enhance be behavior of the interactive mode of VASP code, and leaves the socket-I/O to `VaspInteractive`,
for better maintanance.

## Howto

`patch.py` modifies the `main.F` and `poscar.F` files of the VASP source code.
First, make sure you have uncompressed the VASP source code, modified `makefile` and `makefile.include` accordingly, and VASP program can compile succcessfully. You can test that by running the **Example 1** in the [TLDR](#tl-dr) section.

Locate the directory of the Fortran source codes, such as `~/vasp.6.3.0/src`, and use `patch.py` to apply the changes
```bash
git clone https://github.com/ulissigroup/vasp-interactive.git
cd vasp-interactive/patches
# Change the directory accordingly
python patch.py ~/vasp.6.3.0/src
```

Confirm the new subroutine `INLATT` is patched correctly:
```bash
grep INLATT -5 ~/vasp.6.3.0/src/*.F
```

You should see both `main.F` and `poscar.F` in the grep results.

Then compile VASP using the normal approach:
```bash
cd ~/vasp.6.3.0
make -j8 all
```

## Version compatibility

The interactive mode source code has not been significantly changed in VASP code base since 
version 5.4, meaning the patches is quite likely compatible with even future releases. 
Currently tested versions are:
- 5.4.4pl2
- 6.1.2
- 6.2.x
- 6.3.x

Please contact the maintainer if you have issues applying the patches due to VASP source code change.

## Example builder system

We use the Nvidia HPCSDK docker image as the main builder environment with the following specs:

Container image:
`nvcr.io/nvidia/nvhpc:21.2-devel-cuda_multi-ubuntu20.04`

Additional APT packages
- `make`
- `intel-mkl-full`
- `makedepf90`
- `libfftw3-3`
- `libfftw3-dev`
- `ca-certificates`
- `rsync`
- `unzip`
- `wget`
- `git`
- `python3`





