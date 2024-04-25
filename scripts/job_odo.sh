#!/bin/bash
#SBATCH -A trn023
#SBATCH -J gs-julia-1MPI-1GPU
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 0:10:00
#SBATCH -N 1

date

GS_DIR=/gpfs/wolf2/olcf/trn023/scratch/wgodoy/tutorial/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE settings-files.json

# launch this file with sbatch `$ sbatch job_odo.sh`
