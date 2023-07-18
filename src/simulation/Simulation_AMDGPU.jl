
import AMDGPU

function _init_fields_amdgpu(settings::Settings, mcd::MPICartDomain,
                             T)::Fields{T, 3, <:AMDGPU.ROCArray{T, 3}}
    size_x = mcd.proc_sizes[1]
    size_y = mcd.proc_sizes[2]
    size_z = mcd.proc_sizes[3]

    # should be ones
    u = AMDGPU.ones(T, size_x + 2, size_y + 2, size_z + 2)
    v = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)

    u_temp = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)
    v_temp = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)

    roc_offsets = AMDGPU.ROCArray(mcd.proc_offsets)
    roc_sizes = AMDGPU.ROCArray(mcd.proc_sizes)

    d::Int64 = 6
    minL = Int64(floor(settings.L / 2 - d))
    maxL = Int64(floor(settings.L / 2 + d))

    # @TODO: get ideal grid size and threads
    threads = (16, 16)
    # grid size must be the total number of threads of each direction
    grid = (settings.L, settings.L)

    AMDGPU.wait(AMDGPU.@roc groupsize=threads gridsize=grid _populate_amdgpu!(u,
                                                                              v,
                                                                              roc_offsets,
                                                                              roc_sizes,
                                                                              minL,
                                                                              maxL))

    xy_face_t, xz_face_t, yz_face_t = _get_mpi_faces(size_x, size_y, size_z, T)

    fields = Fields(u, v, u_temp, v_temp, xy_face_t, xz_face_t, yz_face_t)
    return fields
end

function iterate!(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}},
                  settings::Settings,
                  mcd::MPICartDomain) where {T, N}
    _exchange!(fields, mcd)
    # this function is the bottleneck
    _calculate!(fields, settings, mcd)

    # swap the names
    fields.u, fields.u_temp = fields.u_temp, fields.u
    fields.v, fields.v_temp = fields.v_temp, fields.v
end

function _populate_amdgpu!(u, v, offsets, sizes, minL, maxL)
    function is_inside(x, y, z, offsets, sizes)::Bool
        if x < offsets[1] || x >= offsets[1] + sizes[1]
            return false
        end
        if y < offsets[2] || y >= offsets[2] + sizes[2]
            return false
        end
        if z < offsets[3] || z >= offsets[3] + sizes[3]
            return false
        end

        return true
    end

    # local coordinates (this are 1-index already)
    lz = (AMDGPU.workgroupIdx().x - Int32(1)) * AMDGPU.workgroupDim().x +
         AMDGPU.workitemIdx().x
    ly = (AMDGPU.workgroupIdx().y - Int32(1)) * AMDGPU.workgroupDim().y +
         AMDGPU.workitemIdx().y

    if lz <= size(u, 3) && ly <= size(u, 2)

        # get global coordinates
        z = lz + offsets[3] - 1
        y = ly + offsets[2] - 1

        if z >= minL && z <= maxL && y >= minL && y <= maxL
            xoff = offsets[1]

            for x in minL:maxL
                # check if global coordinates for initialization are inside the region
                if !is_inside(x, y, z, offsets, sizes)
                    continue
                end

                # Julia is 1-index, like Fortran :)
                u[x - xoff + 2, ly + 1, lz + 1] = 0.25
                v[x - xoff + 2, ly + 1, lz + 1] = 0.33
            end
        end
    end
end

function _calculate_kernel_amdgpu!(u::AMDGPU.Device.ROCDeviceArray{T, 3, 1},
                                   v::AMDGPU.Device.ROCDeviceArray{T, 3, 1},
                                   u_temp::AMDGPU.Device.ROCDeviceArray{T,
                                                                        3, 1
                                                                        },
                                   v_temp::AMDGPU.Device.ROCDeviceArray{T,
                                                                        3, 1
                                                                        },
                                   sizes::AMDGPU.Device.ROCDeviceArray{Int64, 1
                                                                       },
                                   Du::T, Dv::T, F::T, K::T, noise::T,
                                   dt::T)::Nothing where {T}

    # local coordinates (this are 1-index already)
    k = (AMDGPU.workgroupIdx().x - Int32(1)) * AMDGPU.workgroupDim().x +
        AMDGPU.workitemIdx().x
    j = (AMDGPU.workgroupIdx().y - Int32(1)) * AMDGPU.workgroupDim().y +
        AMDGPU.workitemIdx().y
    i = (AMDGPU.workgroupIdx().z - Int32(1)) * AMDGPU.workgroupDim().z +
        AMDGPU.workitemIdx().z

    # loop through non-ghost cells
    if k == 1 || k >= sizes[3] || j == 1 || j >= sizes[2] || i == 1 ||
       i >= sizes[1]
        return
    end

    @inbounds begin
        u_ijk = u[i, j, k]
        v_ijk = v[i, j, k]

        du = Du * (u[i - 1, j, k] + u[i + 1, j, k] + u[i, j - 1, k] +
              u[i, j + 1, k] + u[i, j, k - 1] + u[i, j, k + 1] -
              6.0 * u_ijk) / 6 - u_ijk * v_ijk^2 +
             F * (1.0 - u_ijk) +
             noise * rand(Distributions.Uniform(-1, 1))
        # + noise * AMDGPU.rand(eltype(u))
        # WIP in AMDGPU.jl, works with CUDA.jl

        dv = Dv *
             (v[i - 1, j, k] + v[i + 1, j, k] + v[i, j - 1, k] +
              v[i, j + 1, k] + v[i, j, k - 1] + v[i, j, k + 1] -
              6.0 * v_ijk) / 6 + u_ijk * v_ijk^2 -
             (F + K) * v_ijk

        # advance the next step
        u_temp[i, j, k] = u_ijk + du * dt
        v_temp[i, j, k] = v_ijk + dv * dt
    end
    return nothing
end

function _calculate!(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}},
                     settings::Settings,
                     mcd::MPICartDomain) where {T, N}
    Du = convert(T, settings.Du)
    Dv = convert(T, settings.Dv)
    F = convert(T, settings.F)
    K = convert(T, settings.k)
    noise = convert(T, settings.noise)
    dt = convert(T, settings.dt)

    roc_sizes = AMDGPU.ROCArray(mcd.proc_sizes)

    threads = (1, 1, 512)
    grid = (settings.L, settings.L, settings.L)

    kernel_time = @elapsed begin

    AMDGPU.wait(AMDGPU.@roc groupsize=threads gridsize=grid _calculate_kernel_amdgpu!(fields.u,
                                                                                      fields.v,
                                                                                      fields.u_temp,
                                                                                      fields.v_temp,
                                                                                      roc_sizes,
                                                                                      Du,
                                                                                      Dv,
                                                                                      F,
                                                                                      K,
                                                                                      noise,
                                                                                      dt)) end

    nx, ny, nz = mcd.proc_sizes[1:3]
    theoretical_fetch_size = 2 *
                             (nx * ny * nz - 8 - 4 * (nx - 2) - 4 * (ny - 2) -
                              4 * (nz - 2)) * sizeof(T)
    theoretical_write_size = 2 * ((nx - 2) * (ny - 2) * (nz - 2)) *
                             sizeof(T)

    println("Theoretical fetch size (GB): ", theoretical_fetch_size * 1e-9)
    println("Theoretical write size (GB):", theoretical_write_size * 1e-9)
    ## Effective memory bandwidth
    datasize = theoretical_fetch_size + theoretical_write_size
    println("Laplacian kernel took: ", kernel_time * 1000,
            " ms effective memory bandwidth: ",
            datasize / kernel_time * 1e-9, " GB/s")
end

function get_fields(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}}) where {T, N}
    u = Array(fields.u)
    u_no_ghost = u[(begin + 1):(end - 1), (begin + 1):(end - 1),
                   (begin + 1):(end - 1)]

    v = Array(fields.v)
    v_no_ghost = v[(begin + 1):(end - 1), (begin + 1):(end - 1),
                   (begin + 1):(end - 1)]
    return u_no_ghost, v_no_ghost
end
