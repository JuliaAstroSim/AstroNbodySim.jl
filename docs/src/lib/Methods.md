# Methods

## Index
```@index
Pages = ["Methods.md"]
```

## Public

### Main functions
```@docs
run
step
preprocessdata
compute_force
compute_potential
output
restart
saverestart
loadrestart
softlen
suggest_softlen
suggest_softlen!
set_softlen!
```

### Particle-Mesh
```@docs
smooth_coef
diff_mat
diff_vec
delta_mat2
delta_mat3
laplace_conv_op
laplace_conv
fft_poisson
fft_poisson!
fdm_poisson
```

### Machine Learning
```@docs
train_cnn_poisson2d
train_cnn_poisson3d
cnn_poisson
```

### Timestep
```@docs
ConstantTimestep
AdaptiveTimestep
init_timesteps
find_next_sync_point_and_drift
advance_and_find_timestep
```

### Energy
```@docs
total_angular_momentum
total_potential
total_momentum
total_kinetic
total_energy
```

### Data processing
```@docs
find_particle
substract_by_id
diff_by_id
```

### MOND
```@docs
nu
nu1
nu2
mond_Milgrom1983
QUMOND_PDM_density
QUMOND_phi
QUMOND_acc
QUMOND
```

### Black Hole
```@docs
r_g
radius_gravity
radius_schwarzschild
pseudoNewtonianPotential
pseudoNewtonianAcc
```

### Elliptic Orbit
```@docs
ellipticSemiMajor
ellipticPeriod
eccentricity
```

### Time scales
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

### Output and Logging
```@docs
outputparallel
setuploggers
LogInfo
DefaultTimer
mkpathIfNotExist
traitstring
```

### Tools
```@docs
write_gadget2_makefile
write_gadget2_param
alter_param
alter_line
add_line
```

## Internal
```@docs

```