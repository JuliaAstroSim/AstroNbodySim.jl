# Methods

```@index
Pages = ["Methods.md"]
```

## Main functions

```@docs
run
step
preprocessdata
mainloop
```

## Force / potential

```@docs
compute_force
compute_force_pairwise
compute_acc
compute_potential
compute_kinetic
softlen
suggest_softlen
suggest_softlen!
set_softlen!
set_softlen_local!
```

## Restart / output

```@docs
output
restart
saverestart
loadrestart
check_restartable
redistribute_data
```

## Particle-Mesh

```@docs
diff_mat
laplace_conv_op
```

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
total_kinetic_local
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
find_particle
substract_by_id
diff_by_id
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
QUMOND
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

```@docs
outputparallel
setuploggers
LogInfo
DefaultTimer
mkpathIfNotExist
traitstring
```

## Utilities

```@docs
write_gadget2_makefile
write_gadget2_param
alter_param
alter_line
add_line
```
