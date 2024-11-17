
# source this file to generate a new IJulia kernel for Perlmutter
module load julia/1.10.4
module load cray-hdf5-parallel
ml use /global/common/software/nersc/julia_hpc_24/modules
module load adios2
export JULIA_ADIOS2_PATH=/global/common/software/nersc/julia_hpc_24/adios2/gnu

# instantiate project packages
julia --project -e 'using Pkg; Pkg.instantiate()'

# capture current environment
julia --project -e 'using IJulia; installkernel("Julia-16-threads", "--project=@.", env=Dict("LD_LIBRARY_PATH" => string(ENV["LD_LIBRARY_PATH"]), "JULIA_NUM_THREADS" => "16", "JULIA_ADIOS2_PATH" => string(ENV["JULIA_ADIOS2_PATH"]) ))'
