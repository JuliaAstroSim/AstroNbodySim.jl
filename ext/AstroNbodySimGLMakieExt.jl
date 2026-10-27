module AstroNbodySimGLMakieExt

# GLMakie 扩展：把实时可视化（Realtime=true）所需的 GLMakie 依赖隔离到
# extension 里。主包不 `using GLMakie`，因此：
#   * 无头 CI（General registry 的 AutoMerge 在无 X11 的 runner 上做
#     `import AstroNbodySim` 检查）不会因 GLFW 需要 DISPLAY 而预编译失败；
#   * 只有显式装了 GLMakie 的用户才会加载它。
#
# 前提：主包的 VisualizationInfo.fig / PlotData 字段是无类型标注的
# （见 src/base/Config.jl 的 mutable struct VisualizationInfo），
# 因此这里不需要主包暴露任何 GLMakie 相关类型。

using AstroNbodySim
using GLMakie
using StructArrays
using Unitful
using PhysicalParticles: pack_xyz, getuLength

# Makie 由 GLMakie 依赖，不在主包 deps 里，故用 GLMakie.Makie 限定访问，
# 避免把 Makie 也加进主包依赖。
const Makie = GLMakie.Makie

function pack_pos(sim::AstroNbodySim.Simulation, ::AstroNbodySim.DirectSum, ::AstroNbodySim.GPU)
    uLength = getuLength(sim.config.units)
    d = StructArray(Array(ustrip.(uLength, sim.simdata.Pos)))
    return pack_xyz(d)
end

function pack_pos(sim::AstroNbodySim.Simulation, ::Union{AstroNbodySim.DirectSum, AstroNbodySim.Tree, AstroNbodySim.FFT, AstroNbodySim.FDM, AstroNbodySim.ML}, ::AstroNbodySim.CPU)
    uLength = getuLength(sim.config.units)
    d = StructArray(ustrip.(uLength, get_all_data(sim).Pos))
    return pack_xyz(d)
end

function estimate_markersize(pos::Array{T,2}) where T
    Max = maximum(pos, dims = 1)
    Min = minimum(pos, dims = 1)
    vol = prod(Max - Min)
    return 0.0005 / cbrt(vol)
end

function makie_scatter(sim::AstroNbodySim.Simulation, GravSolver::AstroNbodySim.Gravity, Device::AstroNbodySim.DeviceType; xlims=(0,1),ylims=(0,1),zlims=(0,1))
    pos = pack_pos(sim, GravSolver, Device)
    sim.visinfo.PlotData = Observable(pos)
    sim.visinfo.fig = GLMakie.scatter(
        sim.visinfo.PlotData;
        markersize = iszero(sim.visinfo.markersize) ? estimate_markersize(pos) : sim.visinfo.markersize,
        markerspace=:data,
        figure = (size = sim.visinfo.size,),
    )

    f,a,p = sim.visinfo.fig
    if a isa Makie.Axis # 2D
        if !isnothing(sim.visinfo.xlims)
            Makie.xlims!(a, sim.visinfo.xlims)
        end
        if !isnothing(sim.visinfo.ylims)
            Makie.ylims!(a, sim.visinfo.ylims)
        end
    elseif a isa Makie.LScene # 3D
        #TODO https://discourse.julialang.org/t/makie-adjusting-axis-scale-limits-in-3d/29231/2
    end
    if !isnothing(sim.visinfo.zlims)
        Makie.zlims!(a, sim.visinfo.zlims)
    end

    display(sim.visinfo.fig)
    sim.visinfo.last_plot_time = time()
end

function update_makie_plot(sim::AstroNbodySim.Simulation, GravSolver::AstroNbodySim.Gravity, Device::AstroNbodySim.DeviceType)
    x, y, z = pack_pos(sim, GravSolver, Device)
    @async sim.visinfo.PlotData[] = pack_pos(sim, GravSolver, Device)
    sim.visinfo.last_plot_time = time()
end

end # module AstroNbodySimGLMakieExt
