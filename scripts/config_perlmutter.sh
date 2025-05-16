#!/bin/bash

# `source GrayScott.jl/scripts/config_perlmutter.sh`

# Packages must be installed in the parallel filesystem 
# and not user's home directory ~/.julia (default) which does not scale
export JULIA_DEPOT_PATH=$SCRATCH/.julia

GS_DIR=$SCRATCH/GrayScott.jl
# remove existing generated Manifest.toml to make this script reusable
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-gnu
# this is just due to an adios2 requirement, parallel HDF5 can be added in the future
module load cray-hdf5-parallel
# module load nvhpc not tested, use cudatoolkit instead
module load cudatoolkit/12.2
ml use /global/common/software/nersc/julia_hpc_24/modules
module load adios2

# module julia 1.11, won't work with julia 1.10
ml use /global/common/software/nersc9/julia/modules/
module load julia/1.11.3

# Required for Julia bindings to point at underlying adios2 modules
export JULIA_ADIOS2_PATH=/global/common/software/nersc/julia_hpc_24/adios2/gnu

# Instantiate the project by installing packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'

# MPIPreferences.jl configures MPI.jl to use the underlying Cray's MPICH
# It will populate LocalPreferences.toml with the correct MPI binary and libraries
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(mpiexec="srun", vendor="cray")'

# Set up CUDA.jl and underlying CUDA runtime
# this will polute Project.toml with CUDA.jl
julia --project=$GS_DIR -e 'using JACC; JACC.set_backend("CUDA")'

# adds an entry to LocalPreferences.toml
julia --project=$GS_DIR -e 'using CUDA; CUDA.set_runtime_version!(v"12.2"; local_toolkit=true)'

# Verify the packages are installed correctly
julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
# JACC.jl and GrayScott.jl won't precompile, but other packages will
julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'
