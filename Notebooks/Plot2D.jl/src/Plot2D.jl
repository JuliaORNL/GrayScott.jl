import ADIOS2
import CairoMakie

# This is NERSC specific path. Please change it to your own path
bp_file = ENV["SCRATCH"]*"/"*ENV["USER"]*"/run001/gs-1MPI-1GPU-64L-F32-JACC-CUDA.bp"

adios = ADIOS2.adios_init_serial()
io = ADIOS2.declare_io(adios, "reader")
reader = ADIOS2.open(io, bp_file, ADIOS2.mode_read)

npixels = 32

# Grid is 64x64x64. Use Tuples () for selection\n",
start = (         16,       16,  32 )
count = (   npixels,  npixels,  1 )

sliceU = Array{Float32, 2}(undef, npixels, npixels)
sliceV = Array{Float32, 2}(undef, npixels, npixels)

steps = ADIOS2.steps(reader)
println("total steps: ", steps)

for step in 1:steps

    ADIOS2.begin_step(reader)

    # These U,V specific lines can be refactored into a function
    varU = ADIOS2.inquire_variable(io, "U")
    @assert varU isa ADIOS2.Variable string("Could not find variable U")
    ADIOS2.set_selection(varU, start, count)

    varV = ADIOS2.inquire_variable(io, "V")
    @assert varU isa ADIOS2.Variable string("Could not find variable V")
    ADIOS2.set_selection(varV, start, count)

    ADIOS2.get(reader, varU, sliceU)
    ADIOS2.get(reader, varV, sliceV)

    # varU, varV are NOT populated (deferred mode)
    ADIOS2.end_step(reader)
    # varU, varV are populated

    # Now plot
    if step % 10 == 0
        println("Showing step", step )
        f = CairoMakie.Figure()
        CairoMakie.heatmap(f[1,1], sliceU)
        CairoMakie.heatmap(f[1,2], sliceV),
        display(f)
        # Save each step figure U,V pair in pdf format
        CairoMakie.save(string("U_V_",step,".pdf"), f)
    end
end

ADIOS2.close(reader)
ADIOS2.adios_finalize(adios)
