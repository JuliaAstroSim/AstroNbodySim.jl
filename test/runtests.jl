#=
cd("AstroNbodySim.jl/test")
=#

using Test

using CSV
using Printf
using DataFrames

using AstroPlot
using Colors
using GLMakie
using CairoMakie
using UnicodePlots
using ColorSchemes

using AstroIO, AstroIC

using Distributed
pids = addprocs(4)

@everywhere using AstroNbodySim

@everywhere using Unitful, UnitfulAstro

@everywhere using PhysicalParticles
@everywhere using PhysicalTrees
@everywhere using PhysicalMeshes

@everywhere astro()

mkpathIfNotExist("output")

using ColorSchemes
colors = ColorSchemes.tab10.colors;

Makie.inline!(true)

using CUDA

include("StaticTest.jl")

include("PM.jl")

IsLocal = false
if IsLocal
    include("Local-Timestep.jl")
    include("Local-Plummer.jl")
    include("Local-ML.jl")
    include("Local-EnergyMomentum.jl")
end