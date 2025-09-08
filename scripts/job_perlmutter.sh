#!/bin/bash
#SBATCH -A trn019
#SBATCH -C gpu
#SBATCH -q shared
#SBATCH -J gs-julia-1MPI-1GPU
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 0:02:00
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --gpus-per-task=1

#export SLURM_CPU_BIND="cores"
GS_DIR=$SCRATCH/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE settings-files.json
