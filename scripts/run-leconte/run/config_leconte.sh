module load "julia/1.9.1"
module load "mpi/openmpi-ppc64le"
module load "gnu/10.2.0"
module load "nvhpc/22.11"

# location of the Project.toml file
export GS_DIR=$HOME/GrayScott.jl

# install packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'
