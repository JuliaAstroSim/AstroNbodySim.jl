# AstroNbodySim.jl

A general-purpose, GPU-accelerated astrophysical N-body simulation library,
released under the GPL-3.0 license.

AstroNbodySim.jl builds on the
broader [JuliaAstroSim](https://github.com/JuliaAstroSim) ecosystem
([PhysicalParticles](https://github.com/JuliaAstroSim/PhysicalParticles.jl),
[PhysicalTrees](https://github.com/JuliaAstroSim/PhysicalTrees.jl),
[PhysicalMeshes](https://github.com/JuliaAstroSim/PhysicalMeshes.jl),
[AstroIC](https://github.com/JuliaAstroSim/AstroIC.jl),
[AstroIO](https://github.com/JuliaAstroSim/AstroIO.jl),
[AstroPlot](https://github.com/JuliaAstroSim/AstroPlot.jl),
[ParallelOperations](https://github.com/JuliaAstroSim/ParallelOperations.jl))
to provide a modular, fully unit-aware, hybrid-parallel simulation toolkit.

```@raw html
<p align="center">
  <img alt="Galactic collision rendered with AstroNbodySim.jl"
       src="https://raw.githubusercontent.com/JuliaAstroSim/AstroNbodySim.jl/main/docs/src/examples/pics/readme/GalacticCollision.gif" />
</p>
```

## Installation

```julia
]add AstroNbodySim
```

or, to install from the latest `main` branch:

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
   help?> run
   search: run

     run(sim::Simulation)

     This function does all the work for you:
     1. `prepare` sim environment and preprocess data
     2. evaluate force, initialize timesteps, output
     3. `mainloop`
   ```

## Gravity solvers at a glance

| Solver | Method | Complexity | Typical use |
| ------ | ------ | ---------- | ----------- |
| [`DirectSum`](@ref) | softened pairwise summation | ``\mathcal{O}(N^2)`` | reference solutions, ``N \lesssim 10^5`` |
| [`Tree`](@ref) | Barnes–Hut octree (Gadget-2 framework) | ``\mathcal{O}(N \log N)`` | large particle counts, ``N \gtrsim 10^5`` |
| [`FFT`](@ref) | spectral Poisson solver on a mesh | ``\mathcal{O}(N_g^3 \log N_g)`` | periodic boxes, 3D galaxy simulations |
| [`FDM`](@ref) | finite-difference Poisson solver | sparse/dense linear solve | 1D/2D verification, non-periodic boundaries |
| [`ML`](@ref) | neural-network Poisson solver (experimental) | network inference | research prototype only |

The particle solvers use the Gadget-2 cubic-spline softening kernel and the
Gadget-2 leapfrog (kick–drift–kick) integrator, with fixed or adaptive
timesteps. See the [Manual](@ref "DirectSum gravity solver") for the formulas,
defaults and parallel execution modes.

## Relation to WaveDM.jl

`AstroNbodySim.jl` also provides the N-body engine used by
[`WaveDM.jl`](https://github.com/JuliaAstroSim/WaveDM.jl). There the
`DirectSum` and `Tree` solvers are used for collisionless baryons (bulge and
stellar disk particles). In the static (`:particles_static`) mode the particle
set is frozen and only supplies a background potential; in the dynamic
(`:particles_dynamic`) mode the particles contribute their potential to the
grid at every SPE step and the wave-DM acceleration is interpolated back to
the particle positions. The `baryon_mode` keyword itself belongs to
`WaveDM.jl`, not to this package. Orbit lookback and tidal-force workflows
built on these solvers are demonstrated in `examples/08-TrajectoryLookback`.

## Supporting and Citing

This software was developed as part of academic research. If you use it in
research, teaching, or other activities, please follow the instructions on the
[CITING page](@ref "Citing AstroNbodySim.jl").

`AstroNbodySim.jl` is the N-body engine of
[`WaveDM.jl`](https://github.com/JuliaAstroSim/WaveDM.jl) and should be cited
together with the following papers:

- **Meng, R. & Dong, X. 2026**, *WaveDM.jl: An Adaptable Simulation Framework for Dynamics of Baryonic and Wave Dark Matter on Galaxy Scales*, Research in Astronomy and Astrophysics — [journal article](http://iopscience.iop.org/article/10.1088/1674-4527/ae9042) · [arXiv:2606.25026](https://arxiv.org/abs/2606.25026) · DOI [10.1088/1674-4527/ae9042](https://doi.org/10.1088/1674-4527/ae9042)

```bibtex
@article{meng2026wavedmjladaptablea,
  title = {{{WaveDM}}.Jl: {{An Adaptable Simulation Framework}} for {{Dynamics}} of {{Baryonic}} and {{Wave Dark Matter}} on {{Galaxy Scales}}},
  shorttitle = {{{WaveDM}}.Jl},
  author = {Meng, Runyu and Dong, Xiaobo},
  year = 2026,
  journal = {Research in Astronomy and Astrophysics},
  issn = {1674-4527},
  doi = {10.1088/1674-4527/ae9042},
  url = {http://iopscience.iop.org/article/10.1088/1674-4527/ae9042},
  eprint = {2606.25026},
  archiveprefix = {arXiv},
  primaryclass = {astro-ph.GA},
  urldate = {2026-07-26}
}
```

Publications built on this code base:

- **Meng, R. & Dong, X. 2026**, *Wave Cold Dark Matter Can Reproduce the Radial Acceleration Relation in X-ray Elliptical Galaxies*, Research in Astronomy and Astrophysics — [journal article](http://iopscience.iop.org/article/10.1088/1674-4527/aeac08) · [arXiv:2610.18969](https://arxiv.org/abs/2610.18969) · DOI [10.1088/1674-4527/aeac08](https://doi.org/10.1088/1674-4527/aeac08)

```bibtex
@article{meng2026wavecold,
  title = {Wave {{Cold Dark Matter Can Reproduce}} the {{Radial Acceleration Relation}} in {{X}}-Ray {{Elliptical Galaxies}}},
  author = {Meng, Runyu and Dong, Xiaobo},
  year = 2026,
  journal = {Research in Astronomy and Astrophysics},
  issn = {1674-4527},
  doi = {10.1088/1674-4527/aeac08},
  url = {http://iopscience.iop.org/article/10.1088/1674-4527/aeac08},
  eprint = {2610.18969},
  archiveprefix = {arXiv},
  primaryclass = {astro-ph.GA},
  urldate = {2026-09-24}
}
```

## Manual Outline

```@contents
Pages = [
    "tutorial/basic.md",
    "tutorial/mond.md",
    "tutorial/psn.md",
    "manual/DirectSum.md",
    "manual/Tree.md",
    "manual/FFT.md",
    "manual/FDM.md",
    "manual/ML.md",
    "lib/Types.md",
    "lib/Methods.md",
]
Depth = 1
```

## Contributing and citing

Issues and pull requests are always welcome — see
[Contributing](@ref "Contributing to AstroNbodySim.jl") for guidelines. If you use
AstroNbodySim.jl in research, please cite the package via
[`CITATION.bib`](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/CITATION.bib);
the same BibTeX entries are also bundled in these docs on the
[Citing page](@ref "Citing AstroNbodySim.jl").
