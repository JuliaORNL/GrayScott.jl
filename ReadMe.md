# GrayScott.jl for Julia for HPC

![GrayScott "U" evolution](https://github.com/JuliaORNL/TutorialJuliaHPC/blob/7a8d4e99a4f4aa05ce74b21f673adee86c9e5512/docs/applications/GrayScott/images/gray-scott.gif)


This is a 3D 7-point stencil code to simulate the following [Gray-Scott
reaction diffusion model](https://doi.org/10.1126/science.261.5118.189) that can run on CPU and GPUs using [JACC.jl](https://github.com/JuliaGPU/JACC.jl):

```
U_t = DU * (U_xx + u_yy + u_zz) - U * V^2 + F * (1 - U) + noise * randn(-1,1)
V_t = DV * (V_xx + v_yy + v_zz) + U * V^2 - (F + k) * V
```

This version contains the following capabilities: 

- CPU/GPU solvers using [JACC.jl](git@github.com:JuliaGPU/JACC.jl.git), GPUs supported: NVIDIA, AMD, Intel and Apple
- Parallel I/O using the [ADIOS2.jl](https://github.com/eschnett/ADIOS2.jl) Julia bindings to [ADIOS2](https://github.com/ornladios/ADIOS2)
- Message passing interface (MPI) using [MPI.jl](https://github.com/JuliaParallel/MPI.jl) Julia bindings to MPI
- Easily switch between float- (Float32) and double- (Float64) precision in the configuration file (Apple GPUs do not support Float64, and Intel GPUs have limited support)
- Data analysis under Notebooks/Plot2D written in Julia and Julia+Jupyter 

## How to run

Pre-requisities:

- Julia version v1.11.0 or greater from [julialang.org/downloads](https://julialang.org/downloads/)

### Run locally 

1. **Set up dependencies**

From the `GrayScott.jl` directory instantiate and use MPI artifact jll (preferred method). 
To use a system provided MPI, see [here](https://juliaparallel.org/MPI.jl/latest/configuration/#using_system_mpi)

```julia

$ julia --project

Julia REPL

julia> ]  

(GrayScott.jl)> instantiate
...
(GrayScott.jl)> <-
julia> using MPIPreferences
julia> MPIPreferences.use_jll_binary()
julia> exit()
```

Julia manages its own packages using [Pkg.jl](https://pkgdocs.julialang.org/v1/), the above would create platform-specific `LocalPreferences.toml` and `Manifest.toml` files.

To use a system provided ADIOS2 library, see [here](https://eschnett.github.io/ADIOS2.jl/dev/#Using-a-custom-or-system-provided-ADIOS2-library). 
It just sets the environment variable `JULIA_ADIOS2_PATH` and build ADIOS2.jl in Julia.


2. **Set up the examples/settings-files.json configuration file**

```
{
    "L": 64,
    "Du": 0.2,
    "Dv": 0.1,
    "F": 0.02,
    "k": 0.048,
    "dt": 1.0,
    "plotgap": 10,
    "steps": 1000,
    "noise": 0.1,
    "output": "gs-julia-1MPI-64L-F32.bp",
    "checkpoint": false,
    "checkpoint_freq": 700,
    "checkpoint_output": "ckpt.bp",
    "restart": false,
    "restart_input": "ckpt.bp",
    "adios_config": "adios2.xml",
    "adios_span": false,
    "adios_memory_selection": false,
    "mesh_type": "image",
    "precision": "Float32",
}
```


3. **Running the simulation**

GrayScott.jl uses JACC.jl for performance portability. Use `JACC.set_backend(backend)`, where `backend = CUDA, AMDGPU, Metal, oneAPI` to setup `LocalPreferences.toml` , see [JACC documentation](https://juliagpu.github.io/JACC.jl/stable/#Supported-backends). To run the simulation:

- `CPU threads`: set CPU threads with `-t` 

    ```
    $ julia --project -t 8 gray-scott.jl examples/settings-files.json
    ```

- `GPU`: 

    ```
    $ julia --project gray-scott.jl examples/settings-files.json
    ```

This would generate an adios2 file from the output entry in the configuration file (e.g. `gs-julia-1MPI-64L-F32.bp`) that can be visualized with ParaView with either the VTX or the FIDES readers.


4. **Running on OLCF Summit and Crusher systems**
The code was tested on the Oak Ridge National Laboratory Leadership Computing Facilities (OLCF): [Summit](https://docs.olcf.ornl.gov/systems/summit_user_guide.html) and [Crusher](https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html). Both are used testing a recent version of Julia [v1.9.0-beta3](https://julialang.org/downloads/#upcoming_release) and a `JULIA_DEPOT_PATH` is required to install packages and artifacts. **DO NOT USE your home directory**. We are providing configuration scripts in `scripts/config_XXX.sh` showing the plumming required for these systems. They need to be executed only once per session from the login nodes. 

To reuse these file the first 3 entries must be modified and run on login-nodes and the PATH poiting at a downloaded Julia binary for the corresponding PowerPC (Summit) and x86-64 (Crusher) architectures. Only "CPU" and "CUDA" backends are supported on Summit, while "CPU" and "AMDGPU" backends are supported on Crusher.

    ```
    # Replace these 3 entries
    PROJ_DIR=/gpfs/alpine/proj-shared/csc383
    export JULIA_DEPOT_PATH=$PROJ_DIR/etc/summit/julia_depot
    GS_DIR=$PROJ_DIR/wgodoy/ADIOS2-Examples/source/julia/GrayScott.jl
    ...
    # and the path 
    export PATH=$PROJ_DIR/opt/summit/julia-1.9.0-beta3/bin:$PATH
    ```

 
## Citation
If you find GrayScott.jl useful, please cite the following [SC'23 WORKS best paper](https://doi.org/10.1145/3624062.3624278).

bib entry:

```
@inproceedings{10.1145/3624062.3624278,
author = {Godoy, William F. and Valero-Lara, Pedro and Anderson, Caira and Lee, Katrina W. and Gainaru, Ana and Ferreira Da Silva, Rafael and Vetter, Jeffrey S.},
title = {Julia as a Unifying End-to-End Workflow Language on the Frontier Exascale System},
year = {2023},
isbn = {9798400707858},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3624062.3624278},
doi = {10.1145/3624062.3624278},
booktitle = {Proceedings of the SC '23 Workshops of The International Conference on High Performance Computing, Network, Storage, and Analysis},
pages = {1989–1999},
numpages = {11},
keywords = {data analysis, exascale, Frontier supercomputer, Jupyter notebooks, end-to-end workflows, High-Performance Computing, Julia, HPC},
location = {Denver, CO, USA},
series = {SC-W '23}
}
```

## Acknowledgements
The work is funded by the US Department of Energy Advanced Scientific Computing Research (ASCR) projects:

- [S4PST]([www.s4pst.org](https://s4pst.org/)) as part of the Next Generation of Scientific Software Technologies (NGSST) ASCR Program.
- NGSST sponsors the Consortium for the Advancement of Scientific Software, [CASS](https://cass.community/)
- ASCR Competitive Portfolios for Advanced Scientific Computing Research, MAGMA/Fairbanks

Former sponsors:

- ASCR Bluestone X-Stack
- The Exascale Computing Project - PROTEAS-TUNE

This research used resources of the Oak Ridge Leadership Computing Facility and the Experimental Computing Laboratory (ExCL) at the Oak Ridge National Laboratory, which is supported by the Office of Science of the U.S. Department of Energy under Contract No. DE-AC05-00OR22725.

Thanks to all the Julia and JuliaGPU community members, packages developers and maintainers for their great work.
