#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 06:00:00



root="/global/u1/t/ttian20/vasp-interactive-test"
export OMP_NUM_THREADS=1
export VASP_COMMAND="srun -n 8 -c 32 --cpu-bind=cores vasp_std"
# vpi is the conda environment prebuilt
module load python
conda activate vpi
pip install -e $root
for ver in  "vasp/5.4.4-cpu" "vasp/6.3.2-cpu"
do
    module load $ver
    echo "Testing VaspInteractive on $ver"
    for f in $root/tests/test*.py
    do
        pytest -svv $f; killall vasp_std || echo ""
    done
    module unload vasp
done
