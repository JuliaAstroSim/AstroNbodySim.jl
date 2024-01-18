module AstroNbodySim


__precompile__(true)

using Reexport
using PrecompileTools

# basic
using Dates
using Printf
using Logging
using Statistics
using LinearAlgebra
using Distributed
using DataFrames
using DocStringExtensions
using Base.Iterators
using Measurements

# user-friendly
using Unitful, UnitfulAstro
using ProgressMeter

# arrays
using StaticArrays
using StructArrays
using SparseArrays
using BangBang

# mesh
using FFTW
using Tullio
using PaddedViews
using OffsetArrays
using IterativeSolvers

# I/O
using FileIO

# Visualize
using GLMakie

# Parallel
using LoopVectorization

# Machine Learning
# using Knet
using Images
using Zygote

# JuliaAstroSim
using ParallelOperations
@reexport using PhysicalParticles
@reexport using PhysicalMeshes
@reexport using PhysicalTrees
@reexport using AstroIO
@reexport using PhysicalFDM
@reexport using AstroSimBase

# GPU
using CUDA
include("base/UtilsGPU.jl")
@hascuda begin
    using CUDA

    println("CUDA-enabled GPU(s) detected:")
    for (gpu, dev) in enumerate(CUDA.devices())
        println(dev)
    end
    CUDA.allowscalar(false)
end

import Base: show, run, step, length, iterate, real, convert
import Unitful: Units, AbstractQuantity

import Distributed: procs
import ParallelOperations: bcast, gather, sum, minimum, maximum
import PhysicalTrees: send_buffer, AbstractTree
import AstroIO: GadgetTypes, HeaderGadget2
import PhysicalMeshes: BoundaryCondition, Vacuum, Periodic, Dirichlet

export
    # Base
    show,
    get_local_data,
    get_all_data,

    # Traits
    AbstractSolverType,
    Gravity, 
        DirectSum, Tree, FDM, FFT, ML,
    Hydro,
        SPH, MHD, FEM, FVM,
    GPUAlgorithm,
        AllPairs, Tiled,
    OutboundLimiter,
        Delete, DS, CoarseMesh,

    # Configs
    SimConfig,
    LogInfo,
    StreamInfo,
    TimeInfo,
    OutputInfo,
    PhysicsInfo,
    VisualizationInfo,
    TreeSimConfig,
    OctreeData,
    Buffer,

    Simulation,

    # Output, log
    output,
    outputparallel,
    setuploggers,
    LogInfo,
    DefaultTimer,
    mkpathIfNotExist,
    traitstring,

    # Parallel
    bcast, procs, gather, sum, minimum, maximum,
    gpuinfo,
    gpu_blocks, gpu_threads,

    # Force
    compute_force,
    softlen,
    suggest_softlen,
    suggest_softlen!, set_softlen!,

    # PM
    smooth_coef,
    diff_mat, diff_vec,
    delta_mat2,
    delta_mat3,
    laplace_conv_op, laplace_conv,
    fft_poisson, fft_poisson!,
    fdm_poisson,

    # ML
    # train_cnn_poisson2d,
    # train_cnn_poisson3d,
    # cnn_poisson,

    # Timestep
    ConstantTimestep, AdaptiveTimestep,
    init_timesteps,
    find_next_sync_point_and_drift,
    advance_and_find_timestep,

    # Energy
    compute_potential,
    total_angular_momentum,
    total_potential,
    total_momentum,
    total_kinetic,
    total_energy,

    # Data Processing
    find_particle,
    substract_by_id,
    diff_by_id,

    # MOND
    nu, nu1, nu2,
    mond_Milgrom1983,
    QUMOND_PDM_density,
    QUMOND_phi,
    QUMOND_acc, QUMOND_acc!,

    # Tools
    write_gadget2_makefile,
    write_gadget2_param,
    alter_param, alter_line,
    add_line,

    ## Black Hole
    r_g, radius_gravity, radius_schwarzschild,
    pseudoNewtonianPotential, pseudoNewtonianAcc,

    ## Elliptic Orbit
    ellipticSemiMajor,
    ellipticPeriod,
    eccentricity,
    
    ## Time scales
    typicalvelocity,
    meandensity,
    crosstime,
    hubbletime,
    relaxtime,
    interactiontime,
    dynamicaltime,
    freefalltime,
    orbitaltime,

    restart, saverestart, loadrestart,
    preprocessdata,
    step,
    run

function __init__()

end

#TODO Internal Type
#Float = Float64 # const?

include("base/Traits.jl")

include("tree/config.jl")
include("base/Config.jl")

include("base/Timestep.jl")

include("base/Parallel.jl")
include("base/Buffer.jl")

include("base/Logging.jl")
include("base/Analyse.jl")
include("base/Timing.jl")

include("base/Output.jl")
include("base/Progress.jl")
include("base/Plot.jl")

include("base/Gravity.jl")
include("base/Potential.jl")
include("base/Energy.jl")

include("base/PostProcessing/Force.jl")
include("base/PostProcessing/Potential.jl")

include("mond/milgrom1983.jl")
include("mond/qumond.jl")

include("tree/gravity.jl")
include("tree/timestep.jl")
include("tree/potential.jl")

include("directsum/gravity.jl")
include("directsum/timestep.jl")
include("directsum/potential.jl")

include("directsumgpu/gravity.jl")
include("directsumgpu/timestep.jl")

include("PM/gravity.jl")
include("PM/fft.jl")
include("PM/cnn.jl")
include("PM/output.jl")
include("PM/timestep.jl")

# include("ML/cnn-poisson/models.jl")
# include("ML/cnn-poisson/dataset.jl")
# include("ML/cnn-poisson/train.jl")

include("restart.jl")
include("run.jl")

include("dataprocessing/find.jl")
include("dataprocessing/diff.jl")

include("tools/TimeScales.jl")
include("tools/SofteningLength.jl")
include("tools/BlackHole.jl")
include("tools/EllipticOrbit.jl")
include("tools/gadget2/makefile.jl")
include("tools/gadget2/param.jl")

include("main.jl")

include("precompile.jl")
end