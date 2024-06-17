#!/bin/bash

# Change the first 3 lines appropriately
PROJ_DIR=/gpfs/wolf2/olcf/trn023/scratch/wgodoy/tutorial
export JULIA_DEPOT_PATH=$PROJ_DIR/.julia
GS_DIR=$PROJ_DIR/GrayScott.jl

# remove existing generated Manifest.toml
rm -f $GS_DIR/Manifest.toml
rm -f $GS_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-gnu-amd/8.4.0
module load cray-mpich
module load adios2/2.9.2
module load julia/1.10.4
# point at current Julia installation until module is available on Odo
# export PATH=/gpfs/wolf2/olcf/trn023/world-shared/opt/julia-1.10.3/bin:$PATH

# Required to point at underlying modules above
export JULIA_ADIOS2_PATH=$OLCF_ADIOS2_ROOT

# MPIPreferences to use Cray's MPICH
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("MPIPreferences")'
julia --project=$GS_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(mpiexec="srun", vendor="cray")'

# Set up underlying rocm
# this will polute Project.toml with AMDGPU.jl
julia --project=$GS_DIR -e 'using Pkg; Pkg.add("AMDGPU")'

export ROCM_PATH=/opt/rocm-5.3.0
# adds an entry to LocalPreferences.toml
julia --project=$GS_DIR -e 'using AMDGPU; AMDGPU.ROCmDiscovery.use_artifacts!(false)'

# Instantiate the project by installing packages in Project.toml
julia --project=$GS_DIR -e 'using Pkg; Pkg.instantiate()'
# Verify the packages are installed correctly
julia --project=$GS_DIR -e 'using Pkg; Pkg.build()'
# JACC.jl and GrayScott.jl won't precompile, but other packages will
julia --project=$GS_DIR -e 'using Pkg; Pkg.precompile()'

# Set jacc-amdgpu as the backend in LocalPreferences.toml
julia --project=$GS_DIR -e 'using GrayScott.GrayScottPreferences; GrayScottPreferences.set_backend("jacc-amdgpu")'
