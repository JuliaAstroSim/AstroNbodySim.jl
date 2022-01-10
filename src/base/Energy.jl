"""
$(TYPEDSIGNATURES)
Compute kinetic energy of particle: `0.5 * p.Mass * p.Vel * p.Vel`
"""
function compute_kinetic(p::AbstractParticle)
    return 0.5 * p.Mass * p.Vel * p.Vel
end

"""
$(TYPEDSIGNATURES)
Sum kinetic energy: `0.5 * data.Mass[i] * data.Vel[i] * data.Vel[i]`
"""
function total_kinetic(data::StructArray)
    s = similar(data.Potential) .* zero(data.Mass[1])
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds s[i] = 0.5 * data.Mass[i] * data.Vel[i] * data.Vel[i]
        end
    end
    return sum(s)
end

"""
$(TYPEDSIGNATURES)
Compute and sum kinetic energy of particles in `data`. Potentials need to be computed in advance.

Return `nothing` if empty.
"""
function total_kinetic_local(sim::Simulation)
    data = get_local_data(sim)
    if iszero(length(data))
        return nothing
    else
        return total_kinetic(data)
    end
end

"""
$(TYPEDSIGNATURES)
Compute kinetic energy of particles on workers and return the sum
"""
function total_kinetic(sim::Simulation)
    return sum(filter(!isnothing, gather(sim, total_kinetic_local)))
end

"""
$(TYPEDSIGNATURES)
Sum potential energy of particles in `data`. Potentials need to be computed in advance.

Return `nothing` if empty.
"""
function total_potential(data::StructArray)
    if iszero(length(data))
        return nothing
    else
        return sum(data.Potential .* data.Mass)
    end
end

total_potential(sim::Simulation) = total_potential(sim, sim.config.solver.grav, sim.config.device.type)

"""
$(TYPEDSIGNATURES)
Total energy (kinetic + potential) of sim data
"""
function total_energy(sim::Simulation)
    k = total_kinetic(sim)
    p = total_potential(sim)
    return k + p
end

"""
$(TYPEDSIGNATURES)
Compute momentum of the system in the direction of `axis`
"""
function total_momentum(sim::Simulation, axis::Symbol)
    data = get_all_data(sim)
    vm = getproperty.(data.Vel, axis) .* data.Mass
    return sum(vm)
end

"""
$(TYPEDSIGNATURES)
Compute total momentum vector of the system
"""
function total_momentum(sim::Simulation)
    data = get_all_data(sim)
    vm = StructArray(data.Vel .* data.Mass)
    return sum(vm)
end

"""
$(TYPEDSIGNATURES)
Compute total angular momentum vector of the system
"""
function total_angular_momentum(sim::Simulation)
    data = get_all_data(sim)
    rvm = data.Mass .* cross.(data.Pos, data.Vel)
    return sum(rvm)
end