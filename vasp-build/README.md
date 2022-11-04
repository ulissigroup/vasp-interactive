# Patch and build VASP binaries with lattice support

This is a brief tutorial to reproduce our patched VASP with an enhanced interactive mode 
that can be used with `VaspInteractive >= v0.1.0`


## **!!!DANGER ZONE!!!**

You are going to change the official VASP source code. Although 
these patches are meant only to modify the behavior of the interactive mode, sudden change / incorrect 
formatting may cause the build to fail or unexpected program behavior. 
Always check the contents of the patches before applying them.


## Do I need these patches?

If you have access to a pre-built VASP which passes our 
[compatibility test](../examples/ex00_vasp_test.py) (see our [report](https://github.com/ulissigroup/vasp-interactive#compatibility-test-fails) for some HPC platforms),
and only need relaxation tasks which do not involve lattice change (e.g. slab geometry optimization, molecular dynamics, etc.), 
`VaspInteractive` should be fully functional.

However, you may consider our patches for the following reasons:

1. No unexpected truncation of output file contents (especially for VASP 5.x)
2. Supporting input of lattice change (e.g. performing lattice optimization, fitting EOS curve)
3. Full compatibility with the iPI socket protocol via `VaspInteractive`

Note we also provide a patch script inspired by the [iPI project](https://github.com/i-pi/i-pi/tree/master/examples/VASP) to add native socket support 
into VASP (5.4 and up). Please check [the advanced topic](#building-native-socket-interface-to-vasp) for more details.




## How-to

### Basic usage
We use the [`compile_vasp.sh`](./compile_vasp.sh) script for the compilation the VASP source code on different linux systems. 
The preprocessing and patching are handled by [`patch.py`](./patch.py). 

We also provide examples of  `makefile` and `makefile.include` under [`makefile_examples`](./makefile_examples/) 
that make us of PGI compilers and Intel MKL libraries. 
Please replace/edit them as needed. 

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

We recommend the container image approach due to better maintenance and reproducibility.
Once you have [installed docker engine](https://docs.docker.com/engine/install/) for your OS, 
you can build the container image and reuse for other platforms with the same CPU architecture.

To use our prepared docker file, simply copy the VASP source code tarball (e.g. `vasp.6.3.0.tgz`)
under the `vasp-build` directory, modify the content of `Dockerfile.example` 
(sections marked with `EDITME`, including VASP version, how to apply the patch, etc.),
 and build the image:

```bash
cd vasp-interactive/vasp-build
# Put vasp.X.Y.Z.tgz under this directory
cp Dockerfile.main.example Dockerfile
# Edit the contents of Dockerfile
docker build . -t <your_image_tag>
```

By default, the compiled VASP binaries are under `/vasp/bin`. You can build production images 
(e.g. including `conda` and `jupyter` environment) on top of the current image. More details see 
the [documentation for docker](https://docs.docker.com/build/)


### Building natively on Linux platforms

`compile_vasp.sh` should be compatible with any Linux system if the correct 
compilation toolkits are chosen.
You can simply adjust the environmental variables passed to `compile_vasp.sh` according to your requirements. Here are a few examples:

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

- **Example 3**: same as Example 2 but with VASP's GPU compatibility via OpenACC

```bash
NCORES=8 ROOT=/tmp VASP_BINARY_PATH=/opt/vasp/bin \
INTERACTIVE_PATCH=patch.py \
./compile_vasp.sh vasp.6.3.0.tgz examples/makefile examples/makefile.include.vasp6-acc
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
can add native support for the socket protocol to VASP. 
We have provided an universal script [`patch_ipi.py`](./patch_ipi.py) for this purpose
(i.e. no need to apply a patch file for each minor version of VASP). 
Since both `patch.py` and `patch_ipi.py` work by matching regex patterns, they can be applied separately to the same source files:
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

### Building intermediate build environment / runtime images (experimental)

For users who wish to directly have a pre-built container image that can directly patch/compile VASP
source codes, we provide two additional docker files, [`Dockerfile.build_env`](./extras/Dockerfile.build_env) 
and [`Dockerfile.runtime`](./extras/Dockerfile.runtime), which are extracted from our main docker file 
to build separate images for building and running VASP, respectively. Build these intermediate images locally first:

```bash
cd vasp-interactive/vasp-build
export DOCKER_BUILDKIT=0
docker build . -f extras/Dockerfile.build_env -t build_env
docker build . --build-arg BUILD_ENV=build_env \
               -f extras/Dockerfile.runtime  -t runtime
```

`build_env` image can be used to produce the VASP binaries on your local file system using the 
[volume mapping](https://docs.docker.com/storage/volumes/). For example, create a container using `build_env`
with an interactive session:

```bash
cd vasp-interactive/vasp-build
docker run -it --rm -v $(pwd):/vasp \
               --workdir /vasp \
               build_env
# Interactive session starts
```
And run the `compile_vasp.sh`:
```bash
INTERACTIVE_PATCH=patch.py ROOT=/tmp VASP_BINARY_PATH=/vasp/bin \
./compile_vasp.sh vasp.X.Y.Z.tgz \
                  makefile_examples/makefile \
                  makefile_examples/makefile.include.vaspX
```
After that, the binaries will be under the `vasp-interactive/vasp-build/bin` directory
of your local file system. 
You can reuse them inside a `runtime` container, for example:
```bash
cd vasp-interactive/vasp-build
docker run --rm -v $(pwd)/bin:/vasp/bin \
           -v <local-vasp-input-files>:/work \
           --workdir /work \
           runtime \
           mpirun -n <cores> --allow-run-as-root /vasp/bin/vasp_std
```
where `local-vasp-input-files` is a directory containing your INCAR, POSCAR etc.

If you belong to the `ulissigroup` Github organization, you will also have access to our pre-built images at:
- `ghcr.io/ulissigroup/vasp-interactive:build_env`
- `ghcr.io/ulissigroup/vasp-interactive:runtime`


Please note due to the EULA licenses of [NVIDIA HPC SDK](https://docs.nvidia.com/hpc-sdk/eula/index.html) 
and [Intel MKL library](https://www.intel.com/content/www/us/en/developer/articles/tool/onemkl-license-faq.html) used,
we are not distributing them publicly. 

 
### Compilation jobs on kubernetes (experimental)

The building process inside the container can also be automated on a kubernetes cluster, as the example
in [`k8s-vasp-build.yaml`](./extras/k8s-vasp-build.yaml) shows:
```bash
# Perform a fresh building
cd vasp-interactive/vasp-build
kubectl delete -f extras/k8s-vasp-build.yaml
kubectl apply -f extras/k8s-vasp-build.yaml
```
You will find the compiled VASP binaries (pristine and patched) under the `workingDir`, which are
the actual ones used for our CI/CD pipelines. Please note the example is only meant to be used for
our internal Kubernetes cluster and should be adapted if you wish to migrate to other systems.







