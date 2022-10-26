# Patch and build VASP binaries with lattice support

This is a brief tutorial to reproduce our patched VASP with enhanced interactive mode 
that can be used with `VaspInteractive >= v0.1.0`


## **!!!DANGER ZONE!!!**

You are going to change the official VASP source code. Although 
these patches are meant only to modify the behavior of the interactive mode, sudden change / incorrect 
formatting may cause build to fail or unexpected program behavior. 
Always check the contents of the patches before applying them.


## Do I need these patches?

If you have access to a pre-built VASP which passes our 
[compatibility test](../examples/ex00_vasp_test.py) (see our [report](https://github.com/ulissigroup/vasp-interactive#compatibility-test-fails) for some HPC platforms),
and only need relaxation tasks which does not involve lattice change (e.g. slab geometry optimization, molecular dynamics, etc.), 
`VaspInteractive` should be fully functional.

However, you may consider apply our patches for the following reasons:

1. No unexpected truncation of output file contents (especially for VASP 5.x)
2. Supporting input of lattice change (e.g. performing lattice optimization, fitting EOS curve)
3. Full compatibility with the iPI socket protocol via `VaspInteractive`

Note we also provide a patch script inspired by the [iPI project](https://github.com/i-pi/i-pi/tree/master/examples/VASP) to add native socket support 
into VASP (5.4 and up). Please check [the advanced topic](#building-native-socket-interface-to-vasp) for more details.

<!-- Note you can also compile the VASP source code to support direct iPI protocol using patches provided with the 
[iPI package](https://github.com/i-pi/i-pi/tree/master/examples/VASP) but currently limited to VASP 5.x. Our patch 
focuses only to enhance be behavior of the interactive mode of VASP code, and leaves the socket-I/O to `VaspInteractive`,
for better maintanance. -->


## How-to

### Basic usage
We use the [`compile_vasp.sh`](./compile_vasp.sh) script for the compilation the VASP source code on different linux systems. 
The preprocessing and patching is handled by [`patch.py`](./patch.py). 

We also provide examples of  `makefile` and `makefile.include` under [`makefile_examples`](./makefile_examples/) 
that make us of PGI compilers and Intel MKL libraries. 
Please replace / edit them as needed. 

As a minimal example, copy your copy of VASP source code `vasp.X.Y.Z.tgz` under 
this directory, and run `compile_vasp.sh`:

```bash
git clone https://github.com/ulissigroup/vasp-interactive.git
cd vasp-interactive/vasp-build
cp <your-local-directory>/vasp.X.Y.Z.tgz ./
# Apply the patch for lattice input and build under /tmp
# your binaries will be under /tmp/vasp.X.Y.Z/bin
INTERACTIVE_PATCH=patch.py ROOT=/tmp \
./compile_vasp.sh vasp.X.Y.Z.tgz \
                  makefile_examples/makefile \
                  makefile_examples/makefile.include.vaspX
```

You may want to check the following sections for more details

### (Recommended) building inside container images

We recommend use the container image approach since you can build the images once and
use for various modern HPC systems (if they support NERSC shifter, Singularity, etc).

There are two pre-built images hosted by this repo:
- `ghcr.io/ulissigroup/vasp-interactive:build_env`
- `ghcr.io/ulissigroup/vasp-interactive:runtime`

`build_env` already contains the `compile_vasp.sh` script which you can build the binaries from a local interactive session:
```bash
cd vasp-interactive/vasp-build
docker run -it --rm -v $(pwd):/work \
               --workdir /work \
               ghcr.io/ulissigroup/vasp-interactive:build_env
# Interactive session starts
```
And run the `compile_vasp.sh`:
```bash
INTERACTIVE_PATCH=patch.py ROOT=/tmp VASP_BINARY_PATH=/work/bin \
./compile_vasp.sh vasp.X.Y.Z.tgz \
                  makefile_examples/makefile \
                  makefile_examples/makefile.include.vaspX
```
The compiled binaries will appear under `vasp-build/bin` of your cloned repo. 
You can reuse them inside the `runtime` container, for example:
```bash
cd vasp-interactive/vasp-build
docker run --rm -v $(pwd)/bin:/vasp \
           -v <local-vasp-input-files>:/work \
           --workdir /work \
           ghcr.io/ulissigroup/vasp-interactive:build_env \
           mpirun -n <cores> --allow-run-as-root /vasp/vasp_std
```
where `local-vasp-input-files` is a directory containing your INCAR, POSCAR etc.

You can also build your own container image based on the `runtime` image. There is a minimal example provided in [`Dockerfile.main.example`](./Dockerfile.main.example). 
Modify and add more components (e.g. `miniconda` packages) as you need and build the image with:
```bash
cd vasp-interactive/vasp-build
cp Dockerfile.main.example Dockerfile
# Edit the contents of Dockerfile
DOCKER_BUILDKIT=0 docker build . -t <your_image_tag>
# Optional: push to registry
docker push <your_image_tag>
```

If you wish to modify and build all images (`build_env`, `runtime`, `your_image_tag`) locally, use the following steps:
```bash
# Edit all Dockerfile.* as needed
export DOCKER_BUILDKIT=0
docker build . -f Dockerfile.build_env -t build_env
docker build . --build-arg BUILD_ENV=build_env \
               -f Dockerfile.runtime  -t runtime
docker build . --build-arg BUILD_ENV=build_env \
               --build-arg RUNTIME=runtime \
               -t <your_image_tag>
# Optional: push to registry
docker push <your_image_tag>
```
which produces local builds of `build_env` and `runtime`.


<!-- The build environment and runtime images for VASP binaries are built from 
[Dockerfile.build_env](Dockerfile.build_env) and [Dockerfile.runtime](Dockerfile.runtime), respectively.
They can be accesed by `ghcr.io/ulissigroup/vasp-interactive:build_env` and  -->


### Building natively on Linux platforms

`compile_vasp.sh` should be compatible with any Linux system provided correct 
compilation toolkits. You can simply adjust the environmental variables according to your requirements. Here are a few examples:

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
For detailed usage please see the source code of `./compile_vasp.sh`. 

The build dependencies vary on different platforms. If you are using Nvidia HPCSDK 
under Ubuntu (e.g. `nvcr.io/nvidia/nvhpc:21.2-devel-cuda_multi-ubuntu20.04`), 
the minimal dependencies are:
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




## Version compatibility

`patch.py` is meant to work for (theoretically) any VASP version > 5.4, 
this is because the interactive mode in VASP source code has not been significantly changed since version 5.4, meaning the patches is quite likely compatible with even future releases.

Please contact us if you have issues applying the patches due to VASP source code change.

## Advanced topics
### Building native socket interface to VASP

The iPI project (https://github.com/i-pi/i-pi) also provides a Fortran plugin that 
can add native support for the socket protocal to VASP. 
We have provided a separate script [`patch_ipi.py`](./patch_ipi.py) to wrap up such process which works indepedent of the VASP version 
(i.e. no need to apply a patch file for each minor version of VASP). 
Since both `patch.py` and `patch_ipi.py` work by matching regex patterns, they can be applied separatedly to the same source files:
```bash
ROOT=/tmp VASP_BINARY_PATH=/work/bin \
INTERACTIVE_PATCH=patch.py \
IPI_PATCH=patch_ipi.py \
./compile_vasp.sh vasp.X.Y.Z.tgz \
                  makefile_examples/makefile \
                  makefile_examples/makefile.include.vaspX
```

The main difference is that
`VaspInteractive` adds the socket communication on top of a local VASP calculator while 
iPI-patch VASP has the socket interface statically compiled. 
In addition, as mentioned in [the introduction](#do-i-need-these-patches),
 `VaspInteractive` provides ad-hoc socket communication compatibility without patching VASP, 
if no lattice change is involved.

### Compilation jobs on kubernetes (internal use)

The building process inside container can also be automated on a kubernetes cluster, as the example
in [`k8s-vasp-build.yaml`](./k8s-vasp-build.yaml) shows:
```bash
# Perform a fresh building
kubectl delete -f k8s-vasp-build.yaml
kubectl apply -f k8s-vasp-build.yaml
```
You will find the compiled VASP binaries (pristine and patched) under the `workingDir`, which are
the actual ones used for our CI/CD pipelines.







