# ML Poisson solver (experimental)

`ML` is a marker type for a convolutional-network Poisson solver. It is
**experimental** and not part of the stable public API: the current source
tree contains an inference hook ([`cnn_poisson`](@ref) in `src/ML/gravity.jl`)
and a training entry point (`cnn_train`) that is still a stub. The network
backend uses `Lux`, which is not declared in the package `Project.toml`, so
the ML path requires a separate Lux installation.

Do not use this solver for production runs; prefer the FFT solver for
periodic boxes and FDM for non-periodic boundary conditions.

## Intended configuration

The constructor passes the model and its parameters through the
`SimConfig` keywords `cnn_model` and `cnn_parameters`:

```julia
sim = Simulation(particles;
                 GravitySolver = ML(),
                 cnn_model     = model,        # Lux model
                 cnn_parameters = parameters,  # trained parameters
                 Nx = Ny = Nz  = 96,
                 xMin = -50u"kpc", xMax = +50u"kpc")
```

Inference is performed by

```julia
cnn_poisson(mesh, unit_projection, solverdata)
```

which evaluates the network on the unit-scaled density field and stores the
result in `mesh.phi`. There is no supported API for training a new network in
this release (`cnn_train` does nothing), and no pretrained checkpoint is
bundled.

## API

```@docs
ML
cnn_poisson
```
