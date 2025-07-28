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
using UnicodePlots
using ColorSchemes

using AstroIO, AstroIC

using Distributed
pids = addprocs(4)

@everywhere using AstroNbodySim
@everywhere using AstroNbodySim.AstroSimBase

@everywhere using Unitful, UnitfulAstro

@everywhere using PhysicalParticles
@everywhere using PhysicalTrees
@everywhere using PhysicalMeshes
@everywhere using AstroNbodySim.PhysicalFDM
@everywhere using AstroNbodySim.PhysicalFFT

@everywhere astro()

mkpathIfNotExist(joinpath(@__DIR__, "output"))

using ColorSchemes
colors = ColorSchemes.tab10.colors;

Makie.inline!(true)

using CUDA


IsLocal = false

include("StaticTest.jl")
include("PM.jl")

if IsLocal
    include("Local-Timestep.jl")
    include("Local-Plummer.jl")
    include("Local-EnergyMomentum.jl")
end
