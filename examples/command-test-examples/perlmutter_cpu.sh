#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -A m2755
#SBATCH -t 00:10:00


export OMP_NUM_THREADS=1

CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
export VASP_COMMAND="srun -u -n 128 --cpu-bind=cores vasp_std"
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi
# conda activate vpi
export PATH=${CONDA_ROOT}/bin:$PATH

# if the $GH_TOKEN variable is correctly set, should return good status
gh auth status

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
CURDIR=`realpath .`
export PYTHONPATH=$CURDIR:$PYTHONPATH
export TEMPDIR=$SCRATCH

for ver in  "vasp/5.4.4-cpu" "vasp/6.3.2-cpu" "vasp-tpc/5.4.4-cpu" "vasp-tpc/6.3.2-cpu"
do
    module load $ver
    echo "Testing VaspInteractive Compatibility on $ver"
    python examples/ex00_vasp_test.py | tee tmp.out
    RES=`sed -n "s/^Test result:\(.*\)$/\1/p" tmp.out`
    echo $ver, $RES >> perlmutter_cpu.txt
    rm tmp.out
    module unload $ver
done

dat=`date +"%Y-%m-%dT%H:%M:%S%z"`
echo "#Last updated: $dat" >> perlmutter_cpu.txt

if [ -z "${GIST_ID}" ]
then
    echo "No gist ID provided, I'll skip the update part"
else
    gh gist edit ${GIST_ID} -a perlmutter_cpu.txt
    echo "Gist file updated!"
fi