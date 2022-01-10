# force from single particle
#=function compute_force_from_particle(pos::AbstractPoint, source::AbstractParticle, G::Number, h::Number)
    dp = source.Pos - pos
    r2 = (dp * dp + h * h) # Softening Kernel
    return G * source.Mass / r2 / sqrt(r2) * dp
end=#

function compute_force_from_particle(pos::AbstractPoint, sourcePos::AbstractPoint, sourceMass::Number, G::Number, h::Number)
    dp = sourcePos - pos
    r2 = dp * dp
    r = sqrt(r2)
    mass = sourceMass
    h_inv = 1.0 / h
    h3_inv = h_inv^3
    
    if r >= h
        fac = mass / (r2 * r)
    else
        u = r * h_inv
        if ustrip(NoUnits, u) < 0.5
            fac = mass * h3_inv * (10.666666666667 + u^2 * (32.0 * u - 38.4))
        else
            fac = mass * h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u^2 - 10.666666666667 * u^3 - 0.066666666667 / u^3)
        end
    end
    
    return dp * fac * (1.0 * G)
end

function compute_force_at_point(particle::AbstractParticle, source::AbstractParticle, G::Number, h::Number)
    return compute_force_from_particle(particle.Pos, source.Pos, source.Mass, G, h)
end

# Compute force from data
#TODO ?Val
#! Do not multi-threading here!
function compute_force_at_point(pos::AbstractPoint, data::StructArray, G::Number, h::Number)
    #acc = [compute_force_from_particle(pos, d, G, h) for d in Iterators.flatten(values(data))]
    #acc = tmap1(d -> compute_force_from_particle(pos, d, G, h), Iterators.flatten(values(data)))
    acc = similar(data.Acc)

    for i in eachindex(data)
        @inbounds acc[i] = compute_force_from_particle(pos, data.Pos[i], data.Mass[i], G, h) 
    end
    return sum(acc)
end
function compute_force_at_point(particle::AbstractParticle, data::StructArray, G::Number, h::Number)
    #acc = [compute_force_from_particle(particle.Pos, d, G, h) for d in Iterators.flatten(values(data))]
    #acc = tmap1(d -> compute_force_from_particle(particle.Pos, d, G, h), Iterators.flatten(values(data)))
    acc = similar(data.Acc)

    for i in eachindex(data)
        @inbounds acc[i] = compute_force_from_particle(particle.Pos, data.Pos[i], data.Mass[i], G, h) 
    end
    return sum(acc)
end
function compute_force_at_point(pos::AbstractPoint, sim::Simulation, h::Number)
    G = sim.config.constants.G
    compute_force_at_point(pos, sim.simdata, G, h)
end
function compute_force_at_point(particle::AbstractParticle, sim::Simulation, h::Number)
    compute_force_at_point(particle.Pos, sim, h)
end


"""
function compute_force_pairwise(p1::AbstractParticle, p2::AbstractParticle, G::Number, h::Number)

    return acc of p1 and p2 in a Tuple
"""
function compute_force_pairwise(p1::AbstractParticle, p2::AbstractParticle, G::Number, h::Number)
    dp = p2.Pos - p1.Pos
    r2 = (dp * dp + h * h)
    acc = G / r2 / sqrt(r2) * dp

    # return p1, p2 acc
    return acc * p2.Mass, -1.0 * acc * p1.Mass
end

function compute_local_force(sim::Simulation, ::DirectSum, ::CPU)
    empty!(sim.buffer.gravToEval)
    empty!(sim.buffer.gravToEvalRecv)

    PosType = typeof(sim.config.ZeroValues.pos)
    OldAccType = typeof(sim.config.ZeroValues.acc.x)
    LenType = typeof(sim.config.ZeroValues.pos.x)

    data = sim.simdata
    G = sim.config.constants.G
    
    #TODO use a macro
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds h = sim.config.grav.ForceSofteningTable[Int(data.Collection[i])]
            @inbounds data.Acc[i] = compute_force_at_point(data.Pos[i], data, G, h)
        end
    end
    

    # Force to compute on remote workers
    grav = Array{GravToEvaluate{PosType, OldAccType, LenType}, 1}()
    for i in eachindex(data)
        push!(grav, GravToEvaluate(i, data.Pos[i], data.Collection[i], data.OldAcc[i], softlen(data.Collection[i], sim)))
    end

    for p in sim.pids
        sim.buffer.gravToEval[p] = grav
    end
end

function compute_pseudo_force(sim::Simulation, ::DirectSum, ::CPU)
    zeroacc = sim.config.ZeroValues.acc
    AccType = typeof(zeroacc)
    for p in sim.pids
        sim.buffer.gravResult[p] = Array{GravResult{AccType}, 1}()
    end
    
    data = sim.simdata
    G = sim.config.constants.G
    for p in keys(sim.buffer.gravToEvalRecv)
        recv = sim.buffer.gravToEvalRecv[p]
        acc = fill(zeroacc, length(recv))

        Threads.@threads for k in 1:Threads.nthreads()
            Head, Tail = split_block(length(recv), k, Threads.nthreads())
            for i in Head:Tail
                @inbounds acc[i] = compute_force_at_point(recv[i].Pos, data, G, recv[i].SoftLen)
            end
        end

        for i in eachindex(recv)
            push!(sim.buffer.gravResult[p], GravResult(recv[i].Index, acc[i], length(sim.simdata)))
        end
    end
end

function combine_force(sim::Simulation, ::DirectSum, ::CPU)
    data = sim.simdata
    for p in keys(sim.buffer.gravResultRecv)
        for recv in sim.buffer.gravResultRecv[p]
            index = recv.Index
            data.Acc[index] += recv.Acc
            data.GravCost[index] += recv.NInteractions
        end
    end
    
    empty!(sim.buffer.gravToEval)
    empty!(sim.buffer.gravToEvalRecv)
    empty!(sim.buffer.gravResult)
    empty!(sim.buffer.gravResultRecv)
end

function compute_force(sim::Simulation, GravSolver::DirectSum, Device::CPU)
    t_FORCE = time_ns()
    bcast(sim, compute_local_force, args = (GravSolver, Device))
    bcast(sim, send_force_buffer)
    bcast(sim, compute_pseudo_force, args = (GravSolver, Device))
    bcast(sim, send_force_result_buffer)
    bcast(sim, combine_force, args = (GravSolver, Device))
    bcast(sim, postprocessing_force, args = (GravSolver, Device))
    add_timer(sim, FORCE, t_FORCE, time_ns())
end

function compute_local_force_at_points(sim::Simulation, SoftLength::Number, GravSolver::DirectSum, Device::CPU)
    pos = sim.buffer.recvbuffer[1]
    sim.buffer.sendbuffer[1] = compute_force_at_point.(pos, sim, SoftLength)
end

function compute_force(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number, GravSolver::DirectSum, Device::CPU) where T<:AbstractPoint3D
    bcast(sim, :buffer, :recvbuffer, Dict(1 => pos))
    bcast(sim, compute_local_force_at_points, args = (SoftLength, GravSolver, Device))

    results = gather(sim, :buffer, :sendbuffer)
    acc = results[1][1]
    for a in results[2:end]
        acc += a[1]
    end

    #TODO MOND - Milgrom 1983
    return acc
end