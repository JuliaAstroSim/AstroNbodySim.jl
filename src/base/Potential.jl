function compute_potential(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number) where T<:AbstractPoint3D
    compute_potential(sim, pos, SoftLength, sim.config.solver.grav, sim.config.device.type)
end

function remove_self_potential_local(sim::Simulation, ::Union{Tree, DirectSum}, ::CPU)
    data = get_local_data(sim)
    G = sim.config.constants.G
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Potential[i] += G * data.Mass[i] / (softlen(data.Collection[i], sim) / 2.8)
        end
    end
end

function remove_self_potential(sim::Simulation, GravSolver::Union{Tree, DirectSum}, Device::CPU)
    bcast(sim, remove_self_potential_local, args = (GravSolver, Device))
end

function compute_potential(sim::Simulation)
    GravSolver = sim.config.solver.grav
    Device = sim.config.device.type
    compute_potential(sim, GravSolver, Device)
    remove_self_potential(sim, GravSolver, Device)
end