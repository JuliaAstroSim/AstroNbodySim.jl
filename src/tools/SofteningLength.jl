
"""
    softlen(c::Collection, table::MVector)
    softlen(p::AbstractParticle, table::MVector)
    softlen(p::AbstractParticle, sim::Simulation)
    softlen(c::Collection, sim::Simulation)

Retrive softening length from table
"""
function softlen(c::Collection, table::MVector)
    return table[Int(c)]
end

function softlen(p::AbstractParticle, table::MVector)
    return softlen(p.Collection, table)
end

function softlen(p::AbstractParticle, sim::Simulation)
    return softlen(p.Collection, sim.config.grav.ForceSofteningTable)
end

function softlen(c::Collection, sim::Simulation)
    return softlen(c, sim.config.grav.ForceSofteningTable)
end

"""
    suggest_softlen(V::Number, N::Int64)
    suggest_softlen(data::Union{Array, StructArray})
    suggest_softlen(data::Union{Array, StructArray}, collection::Collection)
    suggest_softlen(sim::Simulation)

return recommended softening length

See also `suggest_softlen!` and `set_softlen!`
"""
function suggest_softlen(V::Number, N::Int64)
    return cbrt(V / N)
end

function suggest_softlen(data::Union{Array, StructArray})
    e = extent(data)
    return suggest_softlen(PhysicalParticles.volume(e), length(data))
end

function suggest_softlen(data::Union{Array, StructArray}, collection::Collection)
    d = filter(p->p.Collection==collection, data)
    if isempty(d)
        @warn "No $collection particle found."
        return nothing
    else
        return suggest_softlen(d)
    end
end

function suggest_softlen(sim::Simulation)
    return suggest_softlen(get_all_data(sim))
end

"""
    set_softlen_local!(sim::Simulation, c::Collection, h::Number)

Set softening length of collection `c` to `h`.
"""
function set_softlen_local!(sim::Simulation, c::Collection, h::Number)
    sim.config.grav.ForceSofteningTable[Int(c)] = h
    return nothing
end

"""
    set_softlen_local!(sim::Simulation, h::Number)

Set all softening lengths to `h`
"""
function set_softlen_local!(sim::Simulation, h::Number)
    sim.config.grav.ForceSofteningTable .= h
    return nothing
end

function set_softlen!(sim::Simulation, h, ::Union{DirectSum, Tree}, ::CPU)
    bcast(sim, set_softlen_local!; args = (h,))
end

function set_softlen!(sim::Simulation, h, ::Union{FFT, FDM}, ::CPU)
    set_softlen_local!(sim, h)
end

function set_softlen!(sim::Simulation, h, ::DirectSum, ::GPU)
    set_softlen_local!(sim, h)
end

function set_softlen!(sim::Simulation, h)
    set_softlen!(sim, h, sim.config.solver.grav, sim.config.device.type)
end

"""
    suggest_softlen!(sim::Simulation)

Set all softening lengths to suggested value.
"""
function suggest_softlen!(sim::Simulation)
    #TODO: suggest on the basis of different collections of particles
    h = suggest_softlen(sim)
    @info "Setting softening lengths to $h"
    set_softlen!(sim, h)
end