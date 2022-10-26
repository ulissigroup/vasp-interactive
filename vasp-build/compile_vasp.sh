#!/bin/bash
###########################################################################
# A wrapper script for compiling vasp source code
# Usage:
# compile_vasp.sh <vasp.X.Y.Z.tgz> <makefile_path> <makefile_include_path>
# 
# Additional environmental variable:
# ROOT:  root directory to expand the uncompressed file to
# DIRNAME: if not None, replace the original uncompressed name
# NCORES: number of cpus to use. default is 4
# INTERACTIVE_PATCH: if set and not empty, apply the interactive patch
# IPI_PATCH: if set and not empty, apply the ipi patch
# VASP_BINARY_PATH: if set and not empty, sync the executables to this folder
# LIBBEEF: if set and not empty, link libbeef to VASP (e.g. /opt/beef/lib)
##########################################################################

# NVHPC-specific directives
# In some cases the module nvhpc is not evaled at container startup, so add the path directly
NVHPC_MOD="/opt/nvidia/hpc_sdk/modulefiles"
if [[ -d "${NVHPC_MOD}" ]]
then
    export MODULEPATH=${NVHPC_MOD}:$MODULEPATH
    module load nvhpc || true
fi

if [[ "$#" != 3 ]]
then
    echo "Invalid input values."
    echo "Usage"
    echo "compile_vasp.sh <vasp.X.Y.Z.tgz> <makefile_path> <makefile_include_path>"
    exit 1
else
    TGZ=$1
    MKFILE=$2
    MKFILE_INC=$3
fi

if [ -z "$NCORES" ]
then
    NCORES=4
fi
echo "Running with ${NCORES}"

if [ -z "$ROOT" ]
then
    ROOT=`pwd`
else
    ROOT=`realpath $ROOT`
fi


echo $TGZ
echo $MKFILE
echo $MKFILE_INC
echo "Root" $ROOT
echo "INTER" $INTERACTIVE_PATCH
echo "ipi" $IPI_PATCH

echo "Extracting VASP tarball"
# Get the extracted VASP name
# https://unix.stackexchange.com/questions/229504/find-extracted-directory-name-from-tar-file
# This will only show the root path name, not actual uncompressing
XTAR_NAME=$(tar -xvf $TGZ -C $ROOT | head -1 | cut -f1 -d"/")
tar -xf $TGZ -C $ROOT
OLD_ROOT=$ROOT/${XTAR_NAME}
# NEW_ROOT=${OLD_ROOT}${SUFFIX}
if [ -z "$DIRNAME" ]
then
    NEW_ROOT=${OLD_ROOT}
else
    NEW_ROOT=$ROOT/$DIRNAME
    echo ${OLD_ROOT} "-->" ${NEW_ROOT}
    # rm -rf ${NEW_ROOT}
    rsync -ahvP --delete ${OLD_ROOT}/ ${NEW_ROOT}/
fi
echo "Copy makefile and makefile.include"
cp $MKFILE ${NEW_ROOT}/makefile
cp $MKFILE_INC ${NEW_ROOT}/makefile.include


echo "Patch the getshmem.c"
sed -i 's;#include <sys/sem.h>;#include <sys/sem.h>\n#define SHM_NORESERVE 010000;' ${NEW_ROOT}/src/lib/getshmem.c

if [[ ! -z "$LIBBEEF" ]]
then
    echo "Add BEEF-vdW to VASP"
    LIBBEEF=`realpath ${LIBBEEF}`
    if [[ ! -d ${LIBBEEF} ]]
    then
        echo "Directory ${LIBBEEF} does not exist!"
        exit 1
    fi
    sed -i "/^CPP .*/i #Add\nCPP_OPTIONS += -Dlibbeef" ${NEW_ROOT}/makefile.include
	echo -e "\n#LibBEEF support\nLLIBS += -L${LIBBEEF} -lbeef" >> ${NEW_ROOT}/makefile.include
fi

if [[ ! -z "$INTERACTIVE_PATCH" ]]
then
    echo "Patch interactive mode in VASP source code"
    if [[ ! -f "${INTERACTIVE_PATCH}" ]]
    then
        echo "Cannot find the interactive patch script!"
        exit 1
    fi
    python3 ${INTERACTIVE_PATCH} ${NEW_ROOT}/src
fi

if [[ ! -z "${IPI_PATCH}" ]]
then
    echo "Patch iPi interface in VASP source code"
    if [[ ! -f "${IPI_PATCH}" ]]
    then
        echo "Cannot find the iPI patch script!"
        exit 1
    fi
    python3 ${IPI_PATCH} ${NEW_ROOT}/src
fi

sleep 2

echo "Compile VASP"
cd ${NEW_ROOT}
make -j $NCORES all
echo "Compilation finished. Your binaries are at: " ${NEW_ROOT}/bin

# sync binaries
if [[ ! -z "$VASP_BINARY_PATH" ]]
then
    echo "Syncing VASP binaries to " ${VASP_BINARY_PATH}
    rsync -ahvP ${NEW_ROOT}/bin/ ${VASP_BINARY_PATH}/
fi