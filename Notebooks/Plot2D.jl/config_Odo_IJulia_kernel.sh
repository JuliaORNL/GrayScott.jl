#!/bin/bash

# source this file to generate a new IJulia kernel for Odo
# Download your own Julia version 
mkdir opt
cd opt
wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.5-linux-x86_64.tar.gz
tar -xzf julia-1.11.5-linux-x86_64.tar.gz
export PATH=$PWD/julia-1.11.5/bin:$PATH
cd ..

# source this file to generate a new IJulia kernel for Odo
git clone https://github.com/JuliaORNL/GrayScott.jl.git
cd GrayScott.jl/Notebooks/Plot2D.jl

# instantiate project packages
julia --project -e 'using Pkg; Pkg.instantiate()'

# capture current environment
julia --project -e 'using IJulia; installkernel("Julia-16-threads", "--project=@.", env=Dict("LD_LIBRARY_PATH" => string(ENV["LD_LIBRARY_PATH"]), "JULIA_NUM_THREADS" => "16" ))'

cd $HOME
# resulting kernel.json
cat ~/.local/share/jupyter/kernels/Julia-16-threads/kernel.json
