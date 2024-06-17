import ADIOS2
import CairoMakie

bp_file = "/gpfs/wolf2/olcf/trn023/scratch/wgodoy/tutorial/runs/run003/gs-1MPI-1GPU-256L-F32-JACC-AMDGPU.bp"

adios = ADIOS2.adios_init_serial()
io = ADIOS2.declare_io(adios, "reader")

reader = ADIOS2.open(io, bp_file, ADIOS2.mode_read)

# get the number of available steps in the reader, only works for bp files
steps = ADIOS2.steps(reader)

npixels = 256
local_count = 32
# Grid is 1800x1800x1800. Use Tuples () for selection
local_start = convert(Int64, ceil((npixels - local_count) / 2))

start = (local_start, local_start, convert(Int64, ceil(npixels / 2)))
count = (local_count, local_count, 1)

sliceU = Array{Float32, 2}(undef, local_count, local_count)
sliceV = Array{Float32, 2}(undef, local_count, local_count)

for step in 1:steps
    ADIOS2.begin_step(reader)

    varU = ADIOS2.inquire_variable(io, "U")
    @assert varU isa ADIOS2.Variable string("Could not find variable U")
    ADIOS2.set_selection(varU, start, count)

    varV = ADIOS2.inquire_variable(io, "U")
    @assert varV isa ADIOS2.Variable string("Could not find variable V")
    ADIOS2.set_selection(varV, start, count)

    ADIOS2.get(reader, varU, sliceU)
    ADIOS2.get(reader, varV, sliceV)
    ADIOS2.end_step(reader)

    println("Showing step ", step)
    f = CairoMakie.Figure()
    CairoMakie.heatmap(f[1, 1], sliceU)
    CairoMakie.heatmap(f[1, 2], sliceV)
    display(f)
    CairoMakie.save(string("U_V_", step, ".pdf"), f)
end

ADIOS2.close(reader)
ADIOS2.adios_finalize(adios)
