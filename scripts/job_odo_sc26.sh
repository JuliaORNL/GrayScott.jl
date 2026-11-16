#!/bin/bash
#SBATCH -A trn047
#SBATCH -J gs-julia-1MPI-1GPU
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 0:02:00
#SBATCH -N 1

date
#rocminfo

SUB_USER=wgodoy
GS_DIR=/gpfs/wolf2/olcf/trn047/scratch/$USER/sc26/$SUB_USER/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE settings-file-odo.json

# launch this file with sbatch `$ sbatch job_odo.sh`
