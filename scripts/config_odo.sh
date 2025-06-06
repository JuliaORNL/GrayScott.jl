#!/bin/bash

PROJ_DIR=/gpfs/wolf2/olcf/trn036/scratch/$USER
export JULIA_DEPOT_PATH=$PROJ_DIR/.julia
GS_DIR=$PROJ_DIR/GrayScott.jl

# remove existing generated Manifest.toml
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-cray/8.6.0
module load cray-mpich
module load rocm/6.3.1
module load adios2/2.10.0-mpi
module load Core/25.05
module load julia/1.11.0

# AMDGPU.jl requires rocm system libraries (default = /opt/rocm)
export ROCM_PATH=/opt/rocm-6.3.1

# Required to point at underlying modules above
export JULIA_ADIOS2_PATH=$OLCF_ADIOS2_ROOT

# DOWNLOAD JULIA PACKAGES
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'

# Set MPIPreferences to use Cray's MPICH in LocalPreferences.toml
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(mpiexec="srun", vendor="cray")'

# Set JACC AMDGPU back end in LocalPreferences.toml. 
# It will add AMDGPU to Project.toml
julia --project=$GS_DIR -e 'using JACC; JACC.set_backend("AMDGPU")'

# Verify the packages are installed correctly
 julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
 # Precompile dependencies
 julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'
