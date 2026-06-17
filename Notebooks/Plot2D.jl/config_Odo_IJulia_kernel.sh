#!/bin/bash

# Enable Julia v1.11
module load Core/25.05
module load julia/1.11.0

# source this file to generate a new IJulia kernel for Odo
cd GrayScott.jl/Notebooks/Plot2D.jl

# instantiate project packages
julia --project -e 'using Pkg; Pkg.instantiate()'

# capture current environment
julia --project -e 'using IJulia; installkernel("Julia-16-threads", "--project=@.", env=Dict("LD_LIBRARY_PATH" => string(ENV["LD_LIBRARY_PATH"]), "JULIA_NUM_THREADS" => "16" ))'

cd $HOME
# resulting kernel.json
cat ~/.local/share/jupyter/kernels/julia-16-threads-1.11/kernel.json