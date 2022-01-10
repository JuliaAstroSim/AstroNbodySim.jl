function compute_potential(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number) where T<:AbstractPoint3D
    compute_potential(sim, pos, SoftLength, sim.config.solver.grav, sim.config.device.type)
end

function compute_potential(sim::Simulation)
    compute_potential(sim, sim.config.solver.grav, sim.config.device.type)
end