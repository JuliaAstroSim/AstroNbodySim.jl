# Types

```@index
Pages = ["Types.md"]
```

## Core data structures

```@docs
Simulation
SimConfig
LogInfo
StreamInfo
TimeInfo
OutputInfo
PhysicsInfo
VisualizationInfo
GravityConfig
AstroNbodySim.SolverConfig
AstroNbodySim.OutputConfig
AstroNbodySim.TimeConfig
AstroNbodySim.DeviceConfig
```

## Tree and mesh data

```@docs
TreeSimConfig
```

## Trait / marker types

These are singleton types that dispatch the various kernels and
integrators at runtime; the `traitstring` companion returns the
short string used for log files and plot titles.

```@docs
DirectSum
Tree
FDM
FFT
ML
SPH
MHD
FEM
FVM
OutboundLimiter
Delete
DS
CoarseMesh
```
