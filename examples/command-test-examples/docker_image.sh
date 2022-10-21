#!/bin/bash -l
# The command must be run inside one of the ulissigroup images!

SCRATCH=/tmp
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi

if [ -z "$NCORES" ]
then
    NCORES=8
fi
CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
if [ -d "$CONDA_ROOT" ]
then
    export PATH=${CONDA_ROOT}:$PATH
fi

conda activate base
uid=`python -c 'import uuid; print(uuid.uuid4())'`
echo $uid
root=$SCRATCH/vpi-runner/$uid
mkdir -p $root && cd $root
echo "Running tests under $root"
git clone https://github.com/ulissigroup/vasp-interactive.git
cd vasp-interactive
git checkout $GIT_REF
echo "Check to $GIT_REF"
export PYTHONPATH=`realpath .`:$PYTHONPATH
export TEMPDIR=$SCRATCH

res="true"
for ver in "vasp.5.4.4.pl2" "vasp.6.1.2_pgi_mkl" "vasp.6.2.0_pgi_mkl" "vasp.6.3.0_pgi_mkl"
do
    export VASP_COMMAND="mpirun -n $NCORES --bind-to core /opt/$ver/bin/vasp_std"
    echo ${VASP_COMMAND}
    echo "Testing VaspInteractive on $ver"
    # Force using the standard python
    /opt/conda/bin/python examples/ex00_vasp_test.py | tee tmp.out
    RES=`sed -n "s/^Test result:\(.*\)$/\1/p" tmp.out`
    echo $ver, $RES >> ulissi_docker.txt
    rm tmp.out
done

dat=`date +"%Y-%m-%dT%H:%M:%S%z"`
echo "#Last updated: $dat" >> cori_knl.txt

if [ -z "${GIST_ID}" ]
then
    echo "No gist ID provided, I'll skip the update part"
else
    gh gist edit ${GIST_ID} -a ulissi_docker.txt
    echo "Gist file updated!"
fi