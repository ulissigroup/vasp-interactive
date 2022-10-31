#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -A m2755
#SBATCH -t 06:00:00


CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
export TMPDIR=$SCRATCH
export OMP_NUM_THREADS=1
export VASP_COMMAND="srun -u -n 8 -c 32 --cpu-bind=cores vasp_std"
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi
# conda activate $CONDA_ROOT
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
export PYTHONPATH=`realpath .`:$PYTHONPATH
export TEMPDIR=$SCRATCH

res="true"
for ver in "vasp/5.4.4-cpu" "vasp/6.3.2-cpu"
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


gh workflow run perlmutter_cpu_status.yaml -f signal=$res -f jobid=$jobid -f path=$root
