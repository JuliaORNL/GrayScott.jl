#!/bin/bash
#    Begin LSF directives
#BSUB -P trn017
#BSUB -W 00:02
#BSUB -nnodes 1
#BSUB -J gs-julia
#BSUB -o output.%J
#BSUB -e output.%J 
#BSUB -N godoywf@ornl.gov
#    End BSUB directives and begin shell commands


# launch this file with bsub `$ bsub job_ascent.sh` 
# in a directory containing a settings-files.json file
# example file provided in GrayScott.jl/examples/settings-files.json

date
# Modify for your account other than csc383/user
GS_DIR=/gpfs/wolf/proj-shared/trn017/$USER/GrayScott.jl
GS_EXE=$GS_DIR/gray-scott.jl

# GPU runs
jsrun -n 1 -g 1 julia --project=$GS_DIR $GS_EXE settings-files.json

# Multithreaded CPU runs:
# ncores=8
# jsrun -n 1 -c $ncores -g 1 julia --project=$GS_DIR -t $ncores $GS_EXE settings-files.json




