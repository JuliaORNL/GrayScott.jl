#!/bin/bash
#SBATCH -A trn036
#SBATCH -J gs-julia-1MPI-1GPU
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 0:02:00
#SBATCH -N 1

date
#rocminfo

GS_DIR=/gpfs/wolf2/olcf/trn036/scratch/$USER/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE settings-file-odo.json

# launch this file with sbatch `$ sbatch job_odo.sh`
