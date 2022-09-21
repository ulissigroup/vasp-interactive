#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C knl
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 06:00:00

root="/global/u1/t/ttian20/vasp-interactive-test"
export VASP_COMMAND="srun -n 64 -c 4 --cpu-bind=cores vasp_std"
# vpi is the conda environment prebuilt
conda activate vpi
pip install -e $root
for ver in "vasp/5.4.4-knl" "vasp/6.3.0-knl"
do
    module load $ver
    echo "Testing VaspInteractive on $ver"
    for f in $root/tests/test*.py
    do
        pytest -svv $f; killall vasp_std || echo ""
    done
    module unload vasp
done
