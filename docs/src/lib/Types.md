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
SolverConfig
OutputConfig
TimeConfig
DeviceConfig
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
GravityModel
Newton
MOND1983Milgrom
QUMOND
LoggingMode
NormalMode
ProgressMode
SilentMode
TimeIntegration
Euler
Leapfrog
GPUAlgorithm
AllPairs
Tiled
OutboundLimiter
Delete
DS
CoarseMesh
```
