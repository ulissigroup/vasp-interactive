#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 02:00:00

CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi
# conda activate vpi
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
export MODPATH=`realpath .`
# export TMPDIR=$SCRATCH

export SHIFTER_IMAGE=ulissigroup/kubeflow_vasp:beef_vdw

res="true"
for ver in "vasp.5.4.4.pl2" "vasp.6.3.0_pgi_mkl"
do
    echo "Testing VaspInteractive on $ver"
    for f in tests/test*.py
    do
        abs_f=`realpath $f`
        VASP_COMMAND="mpirun -np 32 --bind-to core /opt/${ver}/bin/vasp_std"
        shifter --image=$SHIFTER_IMAGE --env=VASP_COMMAND="$VASP_COMMAND" --env=PYTHONPATH="$MODPATH" --env=TEMPDIR="$SCRATCH" -- pytest -svv ${abs_f}
        # pytest -svv $f
        if [[ $? != 0 ]]
        then
            res="false"
            killall vasp_std
            break
        fi
        killall vasp_std
    done
    module unload vasp
    if [[ $res == "false" ]]
    then
        break
    fi
done

if [[ $res == "true" ]]
then
    echo "All test pass!"
else
    echo "Test fail. See output"
fi


gh workflow run cori_shifter_status.yaml -f signal=$res -f jobid=$jobid -f path=$root

