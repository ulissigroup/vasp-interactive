#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 00:10:00

# For some reason the vasp/6.2.1-hsw and 6.3.2-hsw modules may complain missing dynamic lib
# during runtime. Always use an interactive node and execute this script via bash

CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
# unbuffered output
export VASP_COMMAND="srun -u -n 32 -c 2 --cpu-bind=cores vasp_std"
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi
# #conda activate vpi
export PATH=${CONDA_ROOT}/bin:$PATH

uid=`uuidgen`
root=$SCRATCH/vpi-runner/$uid
mkdir -p $root && cd $root
jobid=${SLURM_JOB_ID}
echo "Job ID $jobid"
echo "Running tests under $root"
gh repo clone $GIT_REPO
cd vasp-interactive
git checkout $GIT_REF
echo "Check to $GIT_REF"
export PYTHONPATH=`realpath .`
export TEMPDIR=$SCRATCH

res="true"
for ver in "vasp/5.4.4-hsw" "vasp/6.2.1-hsw" "vasp/6.3.2-hsw" "vasp-tpc/5.4.4-hsw" "vasp-tpc/6.2.1-hsw" "vasp-tpc/6.3.2-hsw"
# for ver in "vasp/6.3.2-hsw"
do
    module load $ver
    which vasp_std
    echo "Testing VaspInteractive Compatibility on $ver"
    # workdir="/global/u1/t/ttian20/vasp-interactive-test"
    python examples/ex00_vasp_test.py | tee tmp.out
    RES=`sed -n "s/^Test result:\(.*\)$/\1/p" tmp.out`
    echo $ver, $RES >> cori_hsw.txt
    rm tmp.out
    module unload $ver
done

if [ -z "${GIST_ID}" ]
then
    echo "No gist ID provided, I'll skip the update part"
else
    gh gist edit ${GIST_ID} -a cori_hsw.txt
    echo "Gist file updated!"
fi