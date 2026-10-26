# Methods

```@index
Pages = ["Methods.md"]
```

## Main functions

```@docs
run
AstroNbodySim.mainloop
```

## Force / potential

```@docs
compute_force
AstroNbodySim.compute_force_pairwise
compute_acc
softlen
suggest_softlen
suggest_softlen!
set_softlen!
set_softlen_local!
AstroNbodySim.compute_kinetic
```

## Restart / output

`saverestart`, `loadrestart`, `restart`, `check_restartable` and `redistribute_data`
are part of the restart machinery driven by `run` / `Simulation`; they are
internal helpers without user-facing docstrings in this release.

## Particle-Mesh

The Poisson solvers are provided by `PhysicalFFT.jl` and `PhysicalFDM.jl`
(no separate docstrings are attached in this package):

- `fft_poisson(Δ, Len, rho, boundary, device)` and
  `fft_poisson(mesh, G)`.
- `fft_poisson!(mesh, rho, boundary)`.
- `fdm_poisson(mesh, Val(dim), G, mode, device, boundary, sparse)`.

## Machine-Learning solver

The `ML` solver is experimental; see the [ML manual page](@ref "ML Poisson solver (experimental)").
No stable training API is exposed in this release.

## Timestep

```@docs
ConstantTimestep
AdaptiveTimestep
init_timesteps
find_next_sync_point_and_drift
advance_and_find_timestep
```

## Energy & orbits

```@docs
total_angular_momentum
total_potential
total_momentum
total_kinetic
AstroNbodySim.total_kinetic_local
total_energy
ellipticSemiMajor
ellipticPeriod
eccentricity
```

## Time scales

```@docs
typicalvelocity
meandensity
crosstime
hubbletime
relaxtime
interactiontime
dynamicaltime
freefalltime
orbitaltime
```

## Data processing

```@docs
get_all_data
get_local_data
```

## MOND / QUMOND / MOG

```@docs
nu
nu1
nu2
mond_Milgrom1983
QUMOND_PDM_density
QUMOND_phi
QUMOND_acc
QUMOND_acc!
MOG_phi
MOG_acc
MOG_acc!
```

## Pseudo-Newtonian (black-hole) potentials

```@docs
r_g
pseudoNewtonianPotential
pseudoNewtonianAcc
```

## Output & logging

`outputparallel`, `setuploggers`, `DefaultTimer`, `mkpathIfNotExist` and
`traitstring` are internal plumbing around the output/logging pipeline and are
not part of the stable public API in this release (`LogInfo` is documented on
the Types page).

## Utilities

```@docs
write_gadget2_param
alter_param
alter_line
```
