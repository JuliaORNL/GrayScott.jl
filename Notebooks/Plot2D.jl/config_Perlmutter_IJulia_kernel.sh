#!/bin/bash

# source this file to generate a new IJulia kernel for Perlmutter
module load PrgEnv-nvidia
# this is just due to an adios2 requirement, parallel HDF5 can be added in the future
module load cray-hdf5-parallel
# module load nvhpc not tested, use cudatoolkit instead
module load cudatoolkit/12.2
module use /global/common/software/nstaff/blaschke/tutorials/julia-hpc-tutorial-icpp25/nersc/modules
module load adios2
export JULIA_ADIOS2_PATH=/global/common/software/nstaff/blaschke/tutorials/julia-hpc-tutorial-icpp25/nersc/adios2/install/nvidia

module load julia/1.12.1

# instantiate project packages
julia --project -e 'using Pkg; Pkg.instantiate()'

# capture current environment
julia --project -e 'using IJulia; installkernel("Julia-16-threads", "--project=@.", env=Dict("LD_LIBRARY_PATH" => string(ENV["LD_LIBRARY_PATH"]), "JULIA_NUM_THREADS" => "16", "JULIA_ADIOS2_PATH" => string(ENV["JULIA_ADIOS2_PATH"]) ))'
