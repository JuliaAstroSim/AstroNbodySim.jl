# AstroNbodySim.jl

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl)
[![][docs-dev-img]][docs-dev-url]

This is an astrophysical simulation code library under GPL 3.0

[中文Readme](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/master/README中文.md)

## Documentation

- [**Dev**][docs-dev-url] &mdash; *documentation of the in-development version.*

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaastrosim.github.io/AstroNbodySim.jl/dev

For beginners, it is highly recommended to read the [documentation of PhysicalParticles.jl](https://juliaastrosim.github.io/PhysicalParticles.jl/dev/).

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
  - Float16, Float32, Float64, Int128, BigFloat, Measurement, etc.
- Cross-platform: Linux, Windows, MacOS. Easy to deploy
- Hybrid Parallelism: multi-threading, distributed parallelism, GPU acceleration
- Modularity and Versatility: 9 packages, designed for general purposes, highly extentable
- Realtime visualzation (interactive)
- Auto-test workflow

## Features quick view

![Realtime visualization on GPU](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/Plummer.gif)

![Galactic collision](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/GalacticCollision.gif)

![Uncertainty propagation]()

![Autodiff of background potential field]()

![User-difined pipeline: Tidal disruption event (TDE)]()

![User-difined pipeline: TDE accretion history]()



## Supporting and Citing

This software was developed as part of academic research. If you would like to help support it, please star the repository. If you use this software as part of your research, teaching, or other activities, we would be grateful if you could cite the following:

```tex
% arxiv
```

## FAQ

## Package ecosystem

- Basic data structure: [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
- File I/O: [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
- Initial Condition: [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
- Parallelism: [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
- Trees: [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
- Meshes: [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
- Plotting: [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
- Simulation: [AstroNbodySim](https://github.com/JuliaAstroSim/AstroNbodySim.jl)