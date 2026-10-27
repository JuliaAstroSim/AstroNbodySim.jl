# 实时可视化（Realtime=true）相关。
#
# GLMakie 已降为 weakdep + extension（见 ext/AstroNbodySimGLMakieExt.jl），
# 本文件不加载 GLMakie：绘图的具体实现放在 extension 里，这里只做派发。
# 这样无头 CI 上 `import AstroNbodySim` 不会因为 GLFW 需要 X11 display 而预编译失败。

function pack_pos(sim::Simulation, ::DirectSum, ::GPU)
    uLength = getuLength(sim.config.units)
    d = StructArray(Array(ustrip.(uLength, sim.simdata.Pos)))
    return pack_xyz(d)
end

function pack_pos(sim::Simulation, ::Union{DirectSum, Tree, FFT, FDM, ML}, ::CPU)
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

# 取 GLMakie extension；未安装 GLMakie 时返回 nothing。
_makie_ext() = Base.get_extension(@__MODULE__, :AstroNbodySimGLMakieExt)

"""
    makie_scatter(sim, GravSolver, Device; xlims, ylims, zlims)

创建实时可视化窗口。实现在 GLMakie extension 里；若未安装 GLMakie 则报错提示。
"""
function makie_scatter(sim::Simulation, GravSolver::Gravity, Device::DeviceType; xlims=(0,1),ylims=(0,1),zlims=(0,1))
    ext = _makie_ext()
    ext === nothing && error(
        "Realtime visualization requires GLMakie, which is not installed. " *
        "Run `import Pkg; Pkg.add(\"GLMakie\")` (or add it to your project's deps) and try again."
    )
    return ext.makie_scatter(sim, GravSolver, Device; xlims=xlims, ylims=ylims, zlims=zlims)
end

"""
    update_makie_plot(sim, GravSolver, Device)

刷新实时可视化窗口。实现在 GLMakie extension 里；若未安装 GLMakie 则报错提示。
"""
function update_makie_plot(sim::Simulation, GravSolver::Gravity, Device::DeviceType)
    ext = _makie_ext()
    ext === nothing && error(
        "Realtime visualization requires GLMakie, which is not installed. " *
        "Run `import Pkg; Pkg.add(\"GLMakie\")` (or add it to your project's deps) and try again."
    )
    return ext.update_makie_plot(sim, GravSolver, Device)
end
