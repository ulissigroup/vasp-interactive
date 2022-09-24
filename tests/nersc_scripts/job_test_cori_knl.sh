#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C knl
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 02:00:00

CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
export VASP_COMMAND="srun -n 64 -c 4 --cpu-bind=cores vasp_std"
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
export PYTHONPATH=`realpath .`
export TMPDIR=$SCRATCH

res="true"
for ver in "vasp/5.4.4-knl" "vasp/6.3.0-knl"
do
    module load $ver
    echo "Testing VaspInteractive on $ver"
    for f in tests/test*.py
    do
        pytest -svv $f
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


gh workflow run cori_hsw_status.yaml -f signal=$res -f jobid=$jobid -f path=$root

