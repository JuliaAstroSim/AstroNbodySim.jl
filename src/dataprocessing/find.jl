function find_particle(data::StructArray, id::Int)
    for i in eachindex(data)
        if data.ID[i] == id
            return data[i]
        end
    end
    return nothing
end

function find_particle_local(sim::Simulation, id::Int)
    data = get_local_data(sim)
    return find_particle(data, id)
end



function find_particle(sim::Simulation, id::Int)
    #result = gather(sim, find_particle, :data, args = (id))
    for p in sim.pids
        result = sendto(p, find_particle_local, :(registry[$(sim.id)]), AstroNbodySim, args = (id))
        if !isnothing(result)
            return result
        end
    end
    return nothing
end

# TODO function find_particle(sim::Simulation, ids::Vector{Int})