# Trait types
abstract type AbstractSolverType end

abstract type Gravity <: AbstractSolverType end
abstract type Hydro <: AbstractSolverType end

## gravity solvers
"Direct Summation Method"
struct DirectSum <: Gravity end
"Peano-Hilbert Space Filling Octree Method"
struct Tree <: Gravity end
"Finite Differencing Method"
struct FDM <: Gravity end
"Fast Fourier Transform"
struct FFT <: Gravity end
"Machine Learning"
struct ML <: Gravity end

## hydro solvers
"Smoothed Particle Hydrodynamics"
struct SPH <: Hydro end
"Magnetohydrodynamics"
struct MHD <: Hydro end
"Finite Element Method"
struct FEM <: Hydro end
"Finite Volume Method"
struct FVM <: Hydro end

# GPU configuration
abstract type GPUAlgorithm end
struct AllPairs <: GPUAlgorithm end
struct Tiled <: GPUAlgorithm end

# Mesh
"""
    `::OutboundLimiter`. Choose how to handle particles out of the non-periodic simulation box.
    Supported:
        - `Delete`: delete outbound particles
        - `DS`: use direct summation method to compute forces 
        - `CoarseMesh`: Construct a coarse mesh to overlap all particles
"""
abstract type OutboundLimiter end
"Delete outbound particles if they run out of the non-periodic simulation box"
struct Delete <: OutboundLimiter end
"Compute forces using direct summation method if the particles run out of the non-periodic simulation box"
struct DS <: OutboundLimiter end
"Construct a coarse mesh to overlap all particles for non-periodic boundary conditions"
struct CoarseMesh <: OutboundLimiter end
