# AstroNbodySim.jl

A年astrophysical simulation code library under GPL 3.0

## Installation

```julia
]add AstroNbodySim
```
or
```julia
]add https://github.com/JuliaAstroSim/AstroNbodySim.jl
```

**You might need to install [NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-toolkit)**

## Package Features

- Compute with units
- User-friendly
  - Well documented
  - Readable programming
  - Vectorized array operations
  - Dispatch on types for various simulation settings
  - `Float16`, `Float32`, `Float64`, `Int128`, `BigFloat`, `Measurement`, etc.
- Cross-platform: `Linux`, `Windows`, `MacOS`. Easy to deploy
- Hybrid Parallelism: multi-threading, distributed parallelism, GPU acceleration
- Modularity and Versatility: 10+ packages, designed for general purposes, highly extentable
- Realtime visualzation (interactive)
- Auto-test workflow

## Development Guide

1. It is highly recommended to use latest `master` branch of related packages (in JuliaAstroSim), first clone:
   - [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
   - [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
   - [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
   - [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
   - [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
   - [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
   - [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
2. `dev --local [absolute path]` to install packages mentioned above, for example
   ```jl
   pkg> dev --local /home/user/work/AstroNbodySim
   ```
3. `VS Code` as well as its `Julia Extension` is convenient to use, and `Revise.jl` will update changes on the run.

## First Time User Guide

0. Read documentation of `Julia`: https://docs.julialang.org/en/v1/
1. Read documentations and READMEs of related packages
2. Try examples in `AstroNbodySim/examples`. First install packages used in examples by `AstroNbodySim/examples/install_pkgs.jl`.
   The default output directory of all examples is `./output`. The programme would make a directory if the user-defined output path does not exist to avoid unnecessary errors.
3. To interrupt a running simulation, create a file named `stop` in the output directory (same with Gadget2):
   ```sh
   echo > output/stop
   ```
4. It's convenient to check out supported function arguments and keywords by `help?` in `REPL`, for example,
   ```julia
   help?> prepare
   search: prepare preprocessdata

     function prepare(simulation::Simulation)

     Do the following operations:
     1. Say hello
     2. Preprocess data
     3. Check the output directory, make a new one if not exist
     4. Remove "stop" file
     5. Set the global preferred units
     6. Set up logging, timing, profiling and analyzing log files
   ```

## Supporting and Citing

This software was developed as part of academic research. If you would like to help support it, please star the repository. If you use this software as part of your research, teaching, or other activities, we would be grateful if you could cite the following:

```tex
% arxiv
```

## Manual Outline

```@contents
Pages = [
]
Depth = 1
```