# Compute force from sim
function compute_force(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number) where T<:AbstractPoint3D
    compute_force(sim, pos, SoftLength, sim.config.solver.grav, sim.config.device.type)
end

function compute_force(sim::Simulation)
    compute_force(sim, sim.config.solver.grav, sim.config.device.type)
end