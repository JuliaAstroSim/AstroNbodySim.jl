# Tutorial: MOND and QUMOND

`AstroNbodySim.jl` implements Milgrom's modified Newtonian dynamics (MOND) as
a post-processing rescaling of the Newtonian acceleration, and QUMOND as a
mesh-based modified Poisson equation. The model is selected with the
`GravityModel` keyword.

## Selecting a model

| Model | Keyword | Solver support |
| ----- | ------- | -------------- |
| Newtonian | `GravityModel = AstroNbodySim.Newton()` (default) | all solvers |
| Milgrom 1983 | `GravityModel = AstroNbodySim.MOND1983Milgrom()` | `DirectSum`, `Tree`, GPU `DirectSum` |
| QUMOND | `GravityModel = AstroNbodySim.QUMOND()` | `FFT`, `FDM` |

The model types are defined in `AstroSimBase.jl` and are not re-exported by
`AstroNbodySim.jl`; use the qualified names above or `using AstroSimBase`.

The interpolation-function index is controlled by `MOND_nuIndex` (default
`2.0`). The acceleration scale ``a_0`` is not a constructor keyword: it lives
in the simulation constants and is available as
`sim.config.constants.ACC0` (about ``1.2\times10^{-10}\,\mathrm{m\,s^{-2}}``
in SI units).

```julia
sim = Simulation(particles;
                 GravitySolver = DirectSum(),
                 GravityModel  = AstroNbodySim.MOND1983Milgrom(),
                 MOND_nuIndex  = 2.0,
                 ForceSofteningTable = fill(0.1u"kpc", 6))
```

## Milgrom 1983 rescaling

After the Newtonian acceleration ``\mathbf{a}_N`` has been computed, the
particle acceleration is rescaled as

```math
\mathbf{a}_{\mathrm{MOND}}
= \nu\!\left(\frac{\lVert \mathbf{a}_N \rVert}{a_0},\, n\right)
  \mathbf{a}_N ,
```

with

```math
\nu(y,n) = \left[\frac{1}{2}\left(1 + \sqrt{1 + 4\,y^{-n}}\right)\right]^{1/n}.
```

Two standard choices are exposed separately:

- [`nu1`](@ref): ``\nu_1(y) = \tfrac{1}{2}\left(1+\sqrt{1+4/y}\right)``
  (`MOND_nuIndex = 1`).
- [`nu2`](@ref): ``\nu_2(y) = \sqrt{\tfrac{1}{2}\left(1+\sqrt{1+4/y^2}\right)}``
  (`MOND_nuIndex = 2`, the default).

Both satisfy ``\nu \to 1`` for ``y \gg 1``. In the deep-MOND limit
``\nu_2 \sim y^{-1/2}``, so ``a_{\mathrm{MOND}} \sim \sqrt{a_N a_0}``.

The rescaling is applied in `postprocessing_force` for the particle solvers,
using the acceleration of the current step and the integer timeline. It is
not applied to the mesh solvers; use QUMOND there.

## QUMOND

QUMOND removes the non-linearity from the Poisson equation by computing a
"phantom" density from the Newtonian acceleration,

```math
\rho_{\mathrm{PDM}}
= \frac{1}{4\pi G}\,
  \nabla\cdot\left[
    \nu_2\!\left(\frac{\lVert \nabla\phi_N \rVert}{a_0}\right)
    \nabla\phi_N
  \right],
```

and then solving

```math
\nabla^2 \phi = 4\pi G\left(\rho + \rho_{\mathrm{PDM}}\right).
```

In the code the chain is
[`QUMOND_PDM_density`](@ref) → [`QUMOND_phi`](@ref) →
[`QUMOND_acc`](@ref) / [`QUMOND_acc!`](@ref). `QUMOND_acc` returns the
acceleration correction ``-\nabla\phi``; `QUMOND_acc!` adds it in place to
the mesh acceleration. The current implementation solves the QUMOND Poisson
equation with periodic boundaries, so the mesh should be periodic:

```julia
sim_q = Simulation(particles;
                   GravitySolver     = FFT(),
                   GravityModel      = AstroNbodySim.QUMOND(),
                   Nx = Ny = Nz      = 256,
                   xMin = -50u"kpc", xMax = +50u"kpc",
                   BoundaryCondition = Periodic())
run(sim_q)
```

## API

The MOND / QUMOND kernels (`nu`, `nu1`, `nu2`, `mond_Milgrom1983`,
`QUMOND_PDM_density`, `QUMOND_phi`, `QUMOND_acc`, `QUMOND_acc!`, `MOG_phi`,
`MOG_acc`) are documented on the Library → Methods page.
