#=
cd("AstroNbodySim/Benchmark")
=#

using BenchmarkTools
using CairoMakie
using PyPlot
using CUDA
using StructArrays
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using Random
using Serialization
using Test

using LoopVectorization
using Polyester
using ThreadsX
using ThreadTools

using AstroNbodySim
using PhysicalParticles
using PhysicalTrees, PhysicalMeshes
using AstroIO, AstroIC
using AstroPlot
using BenchmarkPlots
using Unitful, UnitfulAstro
using DataFrames

CUDA.allowscalar(false)


if !isdir("output/")
    mkpath("output/")
end

@info "Number of threads: $(Threads.nthreads())"
@info "BLAS threads: $(LinearAlgebra.BLAS.get_num_threads())"

include("DifferencingEquation.jl")
NumData = [2^i for i in 2:12]
benchmark_DifferencingEquation_1D_CPU(NumData; )

NumData = [2^i for i in 2:13]
benchmark_DifferencingEquation_1D_GPU(NumData; )


NumData = [2^i for i in 2:10]
benchmark_DifferencingEquation_1D(NumData)

NumData = [2^i for i in 2:6]
benchmark_DifferencingEquation_2D(NumData)

NumData = collect(4:4:16)
benchmark_DifferencingEquation_3D(NumData)

#=
include("Multi-threading.jl")
=#

include("Scaling.jl")