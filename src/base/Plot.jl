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

function makie_scatter(sim::Simulation, GravSolver::Gravity, Device::DeviceType; xlims=(0,1),ylims=(0,1),zlims=(0,1))
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

function update_makie_plot(sim::Simulation, GravSolver::Gravity, Device::DeviceType)
    x, y, z = pack_pos(sim, GravSolver, Device)
    @async sim.visinfo.PlotData[] = pack_pos(sim, GravSolver, Device)
    sim.visinfo.last_plot_time = time()
end