function DifferencingEquation_CPU_1D(N)
    A = diff_mat(N, 2)
    b = rand(N)
    return A, b
end

function DifferencingEquation_GPU_1D(N)
    A, b = DifferencingEquation_CPU_1D(N)
    return cu(A), cu(b)
end

function DifferencingEquation_CPU_2D(N)
    A = delta_mat2(N,N)
    b = rand(N*N)
    return A, b
end

function DifferencingEquation_GPU_2D(N)
    A, b = DifferencingEquation_CPU_2D(N)
    return cu(A), cu(b)
end

function DifferencingEquation_CPU_3D(N)
    A = delta_mat3(N,N,N)
    b = rand(N*N*N)
    return A, b
end

function DifferencingEquation_GPU_3D(N)
    A, b = DifferencingEquation_CPU_3D(N)
    return cu(A), cu(b)
end

function DifferencingEquation_LinearAlgebra_inv(A, b) #! 1.7 bug: inv of large matrix on CPU -> stackoverflow
    return inv(A) * b
end

function DifferencingEquation_LinearAlgebra_pinv(A, b)
    return pinv(A) * b # BLAS is using all processors by default
end

function DifferencingEquation_LinearAlgebra_left_division(A, b) #! 1.7 bug: inv of large matrix on CPU -> stackoverflow
    return A \ b
end

function DifferencingEquation_LinearAlgebra_left_division_sparse(A, b) #! GPU is not supported
    return sparse(A) \ b
end

# function DifferencingEquation_LinearAlgebra_BLAS(A, b)
#     return A \ b #TODO LinearAlgebra.BLAS.set_num_threads
# end

function DifferencingEquation_turbo(A, b) #! 1.7 bug: inv of large matrix on CPU -> stackoverflow
    return @tturbo A \ b
end


#! IterativeSolvers may not converge!!!

function DifferencingEquation_iterative_cg(A, b)
    return IterativeSolvers.cg(A, b)
end

function DifferencingEquation_iterative_minres(A, b)
    return IterativeSolvers.minres(A, b)
end

function DifferencingEquation_iterative_idrs(A, b) #! Wrong result on GPU
    return IterativeSolvers.idrs(A, b)
end

function DifferencingEquation_iterative_gmres(A, b)  #! Wrong result on CPU. Scalar indexing on GPU
    return IterativeSolvers.gmres(A, b)
end

function DifferencingEquation_iterative_lsmr(A, b)  #! Wrong result on CPU, GPU
    return IterativeSolvers.lsmr(A, b)
end

function DifferencingEquation_iterative_lsqr(A, b)  #! Wrong result on CPU. Scalar indexing on GPU
    return IterativeSolvers.lsqr(A, b)
end

# Stationary methods
function DifferencingEquation_iterative_jacobi(A, b)  #! Wrong result on CPU. Scalar indexing on GPU
    return IterativeSolvers.jacobi(A, b)
end

function DifferencingEquation_iterative_gauss_seidel(A, b)  #! Wrong result on CPU. Scalar indexing on GPU
    return IterativeSolvers.gauss_seidel(A, b)
end



function benchmark_DifferencingEquation_1D_CPU(NumData;
    savefolder = "output",
    title = "Benchmark of solving matrix equations on CPU with $(Threads.nthreads()) threads",
    kw...
)
    fig, df = benchmarkplot(
        [
            DifferencingEquation_LinearAlgebra_pinv,
            DifferencingEquation_LinearAlgebra_left_division_sparse,
            DifferencingEquation_iterative_cg,
            DifferencingEquation_iterative_minres,
            DifferencingEquation_iterative_idrs,
        ],
        [
            "LinearAlgebra_pinv",
            "LinearAlgebra_left_devision_sparse",
            "itreative_cg",
            "itreative_minres",
            "itreative_idrs",
        ],
        DifferencingEquation_CPU_1D,
        NumData;
        savefolder,
        kw...
    )
    Makie.save(joinpath(savefolder, "DifferencingEquation_CPU_1D_$(Threads.nthreads())_threads.png"), fig)
    mv(joinpath(savefolder, "benchmark.csv"), joinpath(savefolder, "DifferencingEquation_CPU_1D_$(Threads.nthreads())_threads.csv"), force = true)
    return fig, df
end

function benchmark_DifferencingEquation_1D_GPU(NumData;
    savefolder = "output",
    title = "Benchmark of solving matrix equations on GPU",
    kw...
)
    fig, df = benchmarkplot(
        [
            #DifferencingEquation_LinearAlgebra_inv, #! Can't do inv on GPU
            DifferencingEquation_LinearAlgebra_pinv,
            DifferencingEquation_LinearAlgebra_left_division,
            DifferencingEquation_turbo,
            DifferencingEquation_iterative_cg,
            DifferencingEquation_iterative_minres,
        ],
        [
            "LinearAlgebra_pinv",
            "LinearAlgebra_left_devision",
            "LoopVectorization",
            "itreative_cg",
            "itreative_minres",
        ],
        DifferencingEquation_GPU_1D,
        NumData;
        savefolder,
        kw...
    )
    Makie.save(joinpath(savefolder, "DifferencingEquation_GPU_1D.png"), fig)
    mv(joinpath(savefolder, "benchmark.csv"), joinpath(savefolder, "DifferencingEquation_GPU_1D.csv"), force = true)
    return fig, df
end

# Plot in one
function benchmark_DifferencingEquation_1D(NumData;
    savefolder = "output",
    title = "Benchmark of solving 1D differencing equations ($(Threads.nthreads()) CPU threads)",
    resolution = (800, 600),
    kw...
)
    Functions = [
        # CPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
        DifferencingEquation_LinearAlgebra_left_division_sparse,

        # GPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
    ]

    Names = [
        #"CPU pinv",
        "CPU left division",
        "CPU left division (sparse)",
        #"GPU pinv",
        "GPU left devision",
    ]

    gen = [
        #DifferencingEquation_CPU_1D,
        DifferencingEquation_CPU_1D,
        DifferencingEquation_CPU_1D,
        #DifferencingEquation_GPU_1D,
        DifferencingEquation_GPU_1D,
    ]

    fig, df = benchmarkplot(
        Functions, Names, gen, NumData;
        savefolder, title, resolution,
        legend = false,
        kw...
    )
    GLMakie.Legend(fig[1,1], fig.scene.children[1].plots[2:end], Names;
        tellheight = false,
        tellwidth = false,
        halign = :left,
        valign = :top,
        margin = (10, 10, 10, 10),
    )
    
    Makie.save(joinpath(savefolder, "BenchmarkDifferencingEquation1D.png"), fig)
    mv(joinpath(savefolder, "benchmark.csv"), joinpath(savefolder, "DifferencingEquation1D.csv"), force = true)
    return fig, df
end

function benchmark_DifferencingEquation_2D(NumData;
    savefolder = "output",
    title = "Benchmark of solving 2D differencing equations ($(Threads.nthreads()) CPU threads)",
    resolution = (800, 600),
    kw...
)
    Functions = [
        # CPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
        DifferencingEquation_LinearAlgebra_left_division_sparse,

        # GPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
    ]

    Names = [
        #"CPU pinv",
        "CPU left division",
        "CPU left division (sparse)",
        #"GPU pinv",
        "GPU left devision",
    ]

    gen = [
        #DifferencingEquation_CPU_2D,
        DifferencingEquation_CPU_2D,
        DifferencingEquation_CPU_2D,
        #DifferencingEquation_GPU_2D,
        DifferencingEquation_GPU_2D,
    ]

    fig, df = benchmarkplot(
        Functions, Names, gen, NumData;
        savefolder, title, resolution,
        legend = false,
        kw...
    )
    GLMakie.Legend(fig[1,1], fig.scene.children[1].plots[2:end], Names;
        tellheight = false,
        tellwidth = false,
        halign = :left,
        valign = :top,
        margin = (10, 10, 10, 10),
    )
    
    Makie.save(joinpath(savefolder, "BenchmarkDifferencingEquation2D.png"), fig)
    mv(joinpath(savefolder, "benchmark.csv"), joinpath(savefolder, "DifferencingEquation2D.csv"), force = true)
    return fig, df
end

function benchmark_DifferencingEquation_3D(NumData;
    savefolder = "output",
    title = "Benchmark of solving 3D differencing equations ($(Threads.nthreads()) CPU threads)",
    resolution = (800, 600),
    kw...
)
    Functions = [
        # CPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
        DifferencingEquation_LinearAlgebra_left_division_sparse,

        # GPU
        #DifferencingEquation_LinearAlgebra_pinv,
        DifferencingEquation_LinearAlgebra_left_division,
    ]

    Names = [
        #"CPU pinv",
        "CPU left division",
        "CPU left division (sparse)",
        #"GPU pinv",
        "GPU left devision",
    ]

    gen = [
        #DifferencingEquation_CPU_3D,
        DifferencingEquation_CPU_3D,
        DifferencingEquation_CPU_3D,
        #DifferencingEquation_GPU_3D,
        DifferencingEquation_GPU_3D,
    ]

    fig, df = benchmarkplot(
        Functions, Names, gen, NumData;
        savefolder, title, resolution,
        legend = false,
        kw...
    )
    GLMakie.Legend(fig[1,1], fig.scene.children[1].plots[2:end], Names;
        tellheight = false,
        tellwidth = false,
        halign = :left,
        valign = :top,
        margin = (10, 10, 10, 10),
    )
    
    Makie.save(joinpath(savefolder, "BenchmarkDifferencingEquation3D.png"), fig)
    mv(joinpath(savefolder, "benchmark.csv"), joinpath(savefolder, "DifferencingEquation3D.csv"), force = true)
    return fig, df
end