# Unit potential from single particle
function compute_unit_potential_from_particle(pos::AbstractPoint, sourcePos::AbstractPoint, sourceMass::Number, G::Number, h::Number)
    dp = sourcePos - pos
    r2 = dp * dp
    r = sqrt(r2)
    mass = sourceMass
    h_inv = 1.0 / h

    if r >= h
        pot = - mass / r * (1.0 * G)
    else
        u = r * h_inv
        if ustrip(NoUnits, u) < 0.5
            wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6))
        else
            wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)))
        end

        pot = mass * h_inv * wp * (1.0 * G)
    end

    if iszero(r)
        pot = pot * 0.0
    end

    return pot
end
function compute_unit_potential_at_point(target::AbstractParticle, source::AbstractParticle, G::Number, h::Number)
    return compute_unit_potential_from_particle(target.Pos, source.Pos, source.Mass, G, h)
end

# Unit potential from data
function compute_unit_potential_at_point(pos::AbstractPoint, data::StructArray, G::Number, h::Number)
    #pot = map(d -> compute_unit_potential_from_particle(pos, d, G, h), Iterators.flatten(values(data)))
    #return sum(pot)
    pot = similar(data.Potential)
    for i in eachindex(data)
        @inbounds pot[i] = compute_unit_potential_from_particle(pos, data.Pos[i], data.Mass[i], G, h)
    end
    return sum(pot)
end
function compute_unit_potential_at_point(particle::AbstractParticle, data::StructArray, G::Number, h::Number)
    #pot = map(d -> compute_unit_potential_from_particle(particle.Pos, d, G, h), Iterators.flatten(values(data)))
    #return sum(pot)
    pot = similar(data.Potential)
    for i in eachindex(data)
        @inbounds pot = compute_unit_potential_from_particle(particle.Pos, data.Pos[i], data.Mass[i], G, h)
    end
    return sum(pot)
end
function compute_unit_potential_at_point(pos::AbstractPoint, sim::Simulation, h::Number)
    G = sim.config.constants.G
    compute_unit_potential_at_point(pos, sim.simdata, G, h)
end
function compute_unit_potential_at_point(particle::AbstractParticle, sim::Simulation, h::Number)
    compute_unit_potential_at_point(particle.Pos, sim, h)
end

function compute_potential_pairwise(p1::AbstractParticle, p2::AbstractParticle, G::Number, h::Number)
    dp = p2.Pos - p1.Pos
    r2 = (dp * dp + h * h)
    pot = - G * p1.Mass * p2.Mass / sqrt(r2)

    return pot
end

function compute_local_potential(sim::Simulation, ::DirectSum, ::CPU)
    empty!(sim.buffer.potToEval)
    empty!(sim.buffer.potToEvalRecv)

    data = sim.simdata
    G = sim.config.constants.G

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Potential[i] = compute_unit_potential_at_point(data.Pos[i], data, G, softlen(data.Collection[i], sim))
        end
    end

    # Force to compute on remote workers
    PosType = typeof(sim.config.ZeroValues.pos)
    pot = Vector{PotentialToEvaluate{PosType}}()
    for i in eachindex(data)
        push!(pot, PotentialToEvaluate(i, data.Pos[i], data.Collection[i], softlen(data.Collection[i], sim)))
    end
    
    for p in sim.pids
        sim.buffer.potToEval[p] = pot
    end
end

function compute_pseudo_potential(sim::Simulation, ::DirectSum, ::CPU)
    empty!(sim.buffer.potToEval)

    potpermass = sim.config.ZeroValues.potpermass
    PotType = typeof(potpermass)
    for p in sim.pids
        sim.buffer.potResult[p] = Vector{PotentialResult{PotType}}()
    end

    data = sim.simdata
    G = sim.config.constants.G
    for p in keys(sim.buffer.potToEvalRecv)
        recv = sim.buffer.potToEvalRecv[p]
        pot = fill(potpermass, length(recv))
        
        Threads.@threads for k in 1:Threads.nthreads()
            Head, Tail = split_block(length(recv), k, Threads.nthreads())
            for i in Head:Tail
                @inbounds pot[i] = compute_unit_potential_at_point(recv[i].Pos, data, G, recv[i].SoftLen)
            end
        end
        
        for i in eachindex(recv)
            push!(sim.buffer.potResult[p], PotentialResult(recv[i].Index, pot[i]))
        end
    end

    empty!(sim.buffer.potToEvalRecv)
end

function combine_potential(sim::Simulation, ::DirectSum, ::CPU)
    empty!(sim.buffer.potResult)
    data = sim.simdata
    
    for p in keys(sim.buffer.potResultRecv)
        for recv in sim.buffer.potResultRecv[p]
            index = recv.Index
            data.Potential[index] += recv.Potential
        end
    end

    # Background potential
    for f in sim.bgpotential
        apply_background_potential(f, data)
    end

    empty!(sim.buffer.potResultRecv)
end

function compute_potential(sim::Simulation, ::DirectSum, ::CPU)
    bcast(sim, compute_local_potential, args = (DirectSum(), CPU()))
    bcast(sim, send_potential_buffer)
    bcast(sim, compute_pseudo_potential, args = (DirectSum(), CPU()))
    bcast(sim, send_potential_result_buffer)
    bcast(sim, combine_potential, args = (DirectSum(), CPU()))
end

function total_potential(sim::Simulation, GravSolver::DirectSum, Device::CPU)
    compute_potential(sim, GravSolver, Device)
    return 0.5 * sum(gather(sim, total_potential, :simdata))
end

function compute_local_potential_at_points(sim::Simulation, SoftLength::Number, ::DirectSum, ::CPU)
    pos = sim.buffer.recvbuffer[1]

    potpermass = sim.config.ZeroValues.potpermass
    # pot = fill(potpermass, size(pos))
    pot = Array{typeof(potpermass), ndims(pos)}(undef, size(pos))
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(pos), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds pot[i] = compute_unit_potential_at_point(pos[i], sim, SoftLength)
        end
    end

    sim.buffer.sendbuffer[1] = pot
end

function compute_potential(sim::Simulation, pos::Array{T,N}, SoftLength::Number, GravSolver::DirectSum, Device::CPU) where T<:AbstractPoint3D where N
    bcast(sim, :buffer, :recvbuffer, Dict(1 => pos))
    bcast(sim, compute_local_potential_at_points, args = (SoftLength, GravSolver, Device))

    results = gather(sim, :buffer, :sendbuffer)
    pot = results[1][1]
    for p in results[2:end]
        pot += p[1]
    end
    return pot
end