#!/bin/bash -l
#SBATCH -A m2755_g
#SBATCH -q debug
#SBATCH -t 0:10:0   
#SBATCH -N 1           
#SBATCH -C gpu
#SBATCH -G 4 

CONDA_ROOT="/global/homes/t/ttian20/.conda/envs/vpi"
export TMPDIR=$SCRATCH
export VASP_COMMAND="srun -n 1 -c 32 --gpu-bind=single:1 -G 4 --cpu-bind=cores vasp_std"
GIT_REPO="ulissigroup/vasp-interactive"
if [ -z "$GIT_REF" ]
then
    GIT_REF="main"
fi
# conda activate $CONDA_ROOT
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
# export PYTHONPATH=`realpath .`
CURDIR=`realpath .`
export PYTHONPATH=$CURDIR:$PYTHONPATH
export TEMPDIR=$SCRATCH

res="true"
for ver in "vasp/6.2.1-gpu" "vasp/6.3.2-gpu" "vasp-tpc/6.2.1-gpu" "vasp-tpc/6.3.2-gpu"
do
    module load $ver
    echo "Testing VaspInteractive on $ver"
    python examples/ex00_vasp_test.py | tee tmp.out
    RES=`sed -n "s/^Test result:\(.*\)$/\1/p" tmp.out`
    echo $ver, $RES >> perlmutter_gpu.txt
    rm tmp.out
done

dat=`date +"%Y-%m-%dT%H:%M:%S%z"`
echo "#Last updated: $dat" >> perlmutter_gpu.txt

if [ -z "${GIST_ID}" ]
then
    echo "No gist ID provided, I'll skip the update part"
else
    gh gist edit ${GIST_ID} -a perlmutter_gpu.txt
    echo "Gist file updated!"
fi