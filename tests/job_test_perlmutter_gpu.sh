#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -q premium
#SBATCH -A m2755
#SBATCH -t 06:00:00




root="/global/u1/t/ttian20/vasp-interactive-test"
export VASP_COMMAND="srun -n1  --cpu-bind=cores --gpu-bind=single:1 -G 4 vasp_std"
# vpi is the conda environment prebuilt
module load python
conda activate vpi
pip install -e $root
for ver in  "vasp/6.2.1-gpu"
do
    module load $ver
    echo "Testing VaspInteractive on $ver"
    for f in $root/tests/test*.py
    do
        pytest -svv $f; killall vasp_std || echo ""
    done
    module unload vasp
done
