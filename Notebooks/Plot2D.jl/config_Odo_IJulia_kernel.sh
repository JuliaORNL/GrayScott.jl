#!/bin/bash

module load julia/1.13.1

export JULIA_DEPOT_PATH=/gpfs/wolf2/olcf/trn047/scratch/$USER/.julia-notebooks

# instantiate project packages
julia --project -e 'using Pkg; Pkg.instantiate()'

# configure IJulia kernel for Odo project
JUPYTER_KERNEL_FILE=$HOME/.local/share/jupyter/kernels/julia-hpc-sc26-1.13/kernel.json

if [ ! -f "$JUPYTER_KERNEL_FILE" ]; then
    echo "Installing IJulia kernel for julia-hpc-sc26"
    julia -e 'using IJulia; installkernel("julia-hpc-sc26", "--project=@.", env=Dict( "USER" => ENV["USER"], "PATH" => ENV["PATH"], "LD_LIBRARY_PATH" => ENV["LD_LIBRARY_PATH"], "JULIA_DEPOT_PATH" => ENV["JULIA_DEPOT_PATH"], "JULIA_NUM_THREADS" => "8"))'
fi

# resulting kernel.json
cat $JUPYTER_KERNEL_FILE
echo "Now log into https://jupyter-open.olcf.ornl.gov/ and select the julia-hpc-sc26 1.13 kernel to run Plot2D.ipynb"
