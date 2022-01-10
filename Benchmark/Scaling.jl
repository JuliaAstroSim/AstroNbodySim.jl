using AstroIC
using DataFrames
using ParameterSpace

function time_sum(Num::Int, solver::String, mode::String)
    OutputDir = "output/$(mode)-$(solver)-$(Num)"
    timings = DataFrame(CSV.File(joinpath(OutputDir, "timing.csv")))
    sum(timings.TOTAL[2:end])
end

function time_minimum(Num::Int, solver::String, mode::String)
    OutputDir = "output/$(mode)-$(solver)-$(Num)"
    timings = DataFrame(CSV.File(joinpath(OutputDir, "timing.csv")))
    minimum(timings.TOTAL[2:end])
end

function time_mean(Num::Int, solver::String, mode::String)
    OutputDir = "output/$(mode)-$(solver)-$(Num)"
    timings = DataFrame(CSV.File(joinpath(OutputDir, "timing.csv")))
    mean(timings.TOTAL[2:end])
end



x = [1, 2, 4, 8]

@warn "Remember to change the path to julia executable"
if Base.Sys.islinux()
    julia_executable = "/mnt/G/linuxwork/julia-1.7.0/usr/bin/julia"
elseif Base.Sys.iswindows()
    julia_executable = "E:\\julia-1.7.0/bin/julia.exe"
end


### Multi-threading
data = generate(PlummerStarCluster(NumSamples = 50000));
write_gadget2("output/scaling.gadget2", data)

df = DataFrame(N = x);
df[:, "DS_multithread_total"]   = zeros(length(x));
df[:, "DS_multithread_min"]     = zeros(length(x));
df[:, "DS_multithread_mean"]    = zeros(length(x));
df[:, "Tree_multithread_total"] = zeros(length(x));
df[:, "Tree_multithread_min"]   = zeros(length(x));
df[:, "Tree_multithread_mean"]  = zeros(length(x));

df[:, "DS_distributed_total"]   = zeros(length(x));
df[:, "DS_distributed_min"]     = zeros(length(x));
df[:, "DS_distributed_mean"]    = zeros(length(x));
df[:, "Tree_distributed_total"] = zeros(length(x));
df[:, "Tree_distributed_min"]   = zeros(length(x));
df[:, "Tree_distributed_mean"]  = zeros(length(x));

# DirectSum multi-threading
for i in x
    run(`$(julia_executable) -t $(i) $(@__DIR__)/scaling/threads-directsum.jl`)
end
df.DS_multithread_total .= time_sum.(x, "DirectSum", "threads")
df.DS_multithread_min .= time_minimum.(x, "DirectSum", "threads")
df.DS_multithread_mean .= time_mean.(x, "DirectSum", "threads")

# Tree multi-threading
for i in x
    run(`$(julia_executable) -t $(i) $(@__DIR__)/scaling/threads-tree.jl`)
end
df.Tree_multithread_total .= time_sum.(x, "Tree", "threads")
df.Tree_multithread_min .= time_minimum.(x, "Tree", "threads")
df.Tree_multithread_mean .= time_mean.(x, "Tree", "threads")

# DirectSum distributed
for i in x
    if i == 1
        run(`$(julia_executable) -t 1 $(@__DIR__)/scaling/distributed-directsum.jl`)
    else
        run(`$(julia_executable) -t 1 -p $(i-1) $(@__DIR__)/scaling/distributed-directsum.jl`)
    end
end
df.DS_distributed_total .= time_sum.(x, "DirectSum", "distributed")
df.DS_distributed_min .= time_minimum.(x, "DirectSum", "distributed")
df.DS_distributed_mean .= time_mean.(x, "DirectSum", "distributed")

# Tree distributed
for i in x
    if i == 1
        run(`$(julia_executable) -t 1 $(@__DIR__)/scaling/distributed-tree.jl`)
    else
        run(`$(julia_executable) -t 1 -p $(i-1) $(@__DIR__)/scaling/distributed-tree.jl`)
    end
end
df.Tree_distributed_total .= time_sum.(x, "Tree", "distributed")
df.Tree_distributed_min .= time_minimum.(x, "Tree", "distributed")
df.Tree_distributed_mean .= time_mean.(x, "Tree", "distributed")

# Save data
CSV.write("output/Scaling.csv", df)

# Plot multi-threading
scene, layout = layoutscene()
axis = layout[1,1] = GLMakie.Axis(
    scene,
    title = "Scaling of mean total time",
    xlabel = "Threads or processes: log2(N)",
    ylabel = "log10(t [ns])",
)
p1 = GLMakie.lines!(axis, log2.(df.N), log10.(df.DS_multithread_mean))
p2 = GLMakie.lines!(axis, log2.(df.N), log10.(df.Tree_multithread_mean))
legend = layout[1,1] = Legend(
    scene,
    [p1,p2],
    [
        "DirectSum multi-threading",
        "Tree multi-threading",
    ];
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    margin = (10, 10, 10, 10),
)
Makie.save("output/ScalingMultiThreading.png", scene)

# Plot distributed
scene, layout = layoutscene()
axis = layout[1,1] = GLMakie.Axis(
    scene,
    title = "Scaling of mean total time",
    xlabel = "Threads or processes: log2(N)",
    ylabel = "log10(t [ns])",
)
p1 = GLMakie.lines!(axis, log2.(df.N), log10.(df.DS_distributed_mean))
p2 = GLMakie.lines!(axis, log2.(df.N), log10.(df.Tree_distributed_mean))
legend = layout[1,1] = Legend(
    scene,
    [p1,p2],
    [
        "DirectSum distributed",
        "Tree distributed",
    ];
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    margin = (10, 10, 10, 10),
)
Makie.save("output/ScalingDistributed.png", scene)

# Plot multi-threading and distributed
scene, layout = layoutscene(;resolution=(800,450))
axis = layout[1,1] = GLMakie.Axis(
    scene,
    title = "Scaling of mean total time",
    xlabel = "Threads or processes: log2(N)",
    ylabel = "log10(t [ns])",
)
p1 = GLMakie.lines!(axis, log2.(df.N), log10.(df.DS_multithread_mean))
p2 = GLMakie.lines!(axis, log2.(df.N), log10.(df.DS_distributed_mean))
p3 = GLMakie.lines!(axis, log2.(df.N), log10.(df.Tree_multithread_mean))
p4 = GLMakie.lines!(axis, log2.(df.N), log10.(df.Tree_distributed_mean))
legend = layout[1,1] = Legend(
    scene,
    [p1,p2,p3,p4],
    [
        "DirectSum multi-threading",
        "DirectSum distributed",
        "Tree multi-threading",
        "Tree distributed",
    ];
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    margin = (10, 10, 10, 10),
)
Makie.save("output/Scaling.png", scene)


### Number of particles
#NumData = [100, 1000, 10000, 100000, 200000, 300000]
NumData = [100, 1000, 10000, 100000]
df = DataFrame(N = NumData);
df[:, "DS_total"]   = zeros(length(NumData));
df[:, "DS_min"]     = zeros(length(NumData));
df[:, "DS_mean"]    = zeros(length(NumData));
df[:, "DSGPU_total"]   = zeros(length(NumData));
df[:, "DSGPU_min"]     = zeros(length(NumData));
df[:, "DSGPU_mean"]    = zeros(length(NumData));
df[:, "Tree_total"] = zeros(length(NumData));
df[:, "Tree_min"]   = zeros(length(NumData));
df[:, "Tree_mean"]  = zeros(length(NumData));


run(`$(julia_executable) -t 8 $(@__DIR__)/scaling/num-gpu-directsum.jl $(NumData)`)
run(`$(julia_executable) -t 8 $(@__DIR__)/scaling/num-threads-directsum.jl $(NumData)`)
run(`$(julia_executable) -t 8 $(@__DIR__)/scaling/num-threads-tree.jl $(NumData)`)


df.DS_total .= time_sum.(NumData, "DirectSum", "threads")
df.DS_min .= time_minimum.(NumData, "DirectSum", "threads")
df.DS_mean .= time_mean.(NumData, "DirectSum", "threads")
df.DSGPU_total .= time_sum.(NumData, "DirectSum", "GPU")
df.DSGPU_min .= time_minimum.(NumData, "DirectSum", "GPU")
df.DSGPU_mean .= time_mean.(NumData, "DirectSum", "GPU")
df.Tree_total .= time_sum.(NumData, "Tree", "threads")
df.Tree_min .= time_minimum.(NumData, "Tree", "threads")
df.Tree_mean .= time_mean.(NumData, "Tree", "threads")
CSV.write("output/ScalingNumber.csv", df)

#= tuning plots
df = DataFrame(CSV.File("output/ScalingNumber.csv"))
=#

scene, layout = layoutscene(; resolution=(800,450))
axis = layout[1,1] = GLMakie.Axis(
    scene,
    title = "Scaling of mean total time",
    xlabel = "log10(N)",
    ylabel = "log10(t [ns])",
)
p1 = GLMakie.lines!(axis, log10.(df.N), log10.(df.DS_mean))
p2 = GLMakie.lines!(axis, log10.(df.N), log10.(df.DSGPU_mean))
p3 = GLMakie.lines!(axis, log10.(df.N), log10.(df.Tree_mean))
legend = layout[1,1] = Legend(
    scene,
    [
        p1,
        p2,
        p3,
    ],
    [
        "DirectSum multi-threading",
        "DirectSum GPU",
        "Tree multi-threading",
    ];
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
    margin = (10, 10, 10, 10),
)
Makie.save("output/ScalingNumber.png", scene)