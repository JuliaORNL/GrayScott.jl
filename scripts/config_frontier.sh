#!/bin/bash

# Change the first 4 lines appropriately
PROJ_DIR=/lustre/orion/proj-shared/csc383/wgodoy
export JULIA_DEPOT_PATH=$PROJ_DIR/../julia_depot
GS_DIR=$PROJ_DIR/GrayScott.jl
JULIA_DIR=/lustre/orion/proj-shared/csc383/opt/julia-1.9.2

# remove existing generated Manifest.toml
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-cray/8.3.3 # has required gcc
module load cray-mpich
module load rocm/5.4.0
module load adios2/2.8.3 # only works with PrgEnv-cray
# module load julia/1.9.0 # currently not available
# patch until julia/1.9.0 is available
export PATH=$JULIA_DIR/bin:$PATH

# Required to point at underlying modules above
export JULIA_AMDGPU_DISABLE_ARTIFACTS=1
# Required to enable underlying ADIOS2 library from loaded module
export JULIA_ADIOS2_PATH=$OLCF_ADIOS2_ROOT

# MPIPreferences to use Cray's MPICH
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("MPIPreferences")'
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["libmpi_cray"], mpiexec="srun")'

# Adds a custom branch in case the development version is needed (for devs to test new features)
julia --project=$GS_DIR -e 'using Pkg; Pkg.add(url="https://github.com/utkarsh530/AMDGPU.jl.git", rev="u/random")'

# Instantiate the project by installing packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'

# Verify the packages are installed correctly
julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'
