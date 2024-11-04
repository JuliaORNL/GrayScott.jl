#!/bin/bash

# Using a predefined $SCRATCH/$USER directory
GS_DIR=$SCRATCH/$USER/GrayScott.jl
# remove existing generated Manifest.toml
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-gnu
module load cray-hdf5-parallel
#module load nvhpc
module load cudatoolkit/12.2
ml use /global/common/software/nersc/julia_hpc_24/modules
module load adios2
module load julia/1.10.4

# Required to point at underlying modules above
export JULIA_ADIOS2_PATH=/global/common/software/nersc/julia_hpc_24/adios2/gnu

# MPIPreferences to use Cray's MPICH
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("MPIPreferences")'
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(mpiexec="srun", vendor="cray")'

# Set up CUDA.jl and underlying CUDA runtime
# this will polute Project.toml with CUDA.jl
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("CUDA")'

# adds an entry to LocalPreferences.toml
julia --project=$GS_DIR -e 'using CUDA; CUDA.set_runtime_version!(v"12.2"; local_toolkit=true)'

# Instantiate the project by installing packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'
# Verify the packages are installed correctly
julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
# JACC.jl and GrayScott.jl won't precompile, but other packages will
julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'

# Set jacc-cuda as the backend in LocalPreferences.toml
# To change back end to CPU modify LocalPreferences.toml or use these commands
julia --project=$GS_DIR -e 'using GrayScott.GrayScottPreferences; GrayScottPreferences.set_backend("jacc-cuda")'
