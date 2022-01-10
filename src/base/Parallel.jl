const registry=Dict{Pair{Int64, Int64},Any}()

let DID::Int = 1
    global next_simulationid
    next_simulationid() = (id = DID; DID += 1; Pair(myid(), id))
end

"""
    next_simulationid()

Produces an incrementing ID that will be used for simulations.
"""
next_simulationid # Document

procs(sim::Simulation) = sim.pids

bcast(sim::Simulation, expr, data, mod::Module = AstroNbodySim) = bcast(sim.pids, :(registry[$(sim.id)].$expr), data, mod)
bcast(sim::Simulation, expr1, expr2, data, mod::Module = AstroNbodySim) = bcast(sim.pids, :(registry[$(sim.id)].$expr1.$expr2), data, mod)
bcast(sim::Simulation, f::Function, expr, mod::Module = AstroNbodySim; args...) = bcast(sim.pids, f, :(registry[$(sim.id)].$expr), mod; args...)
bcast(sim::Simulation, f::Function, expr1, expr2, mod::Module = AstroNbodySim; args...) = bcast(sim.pids, f, :(registry[$(sim.id)].$expr1.$expr2), mod; args...)
bcast(sim::Simulation, f::Function, mod::Module = AstroNbodySim; args...) = bcast(sim.pids, f, :(registry[$(sim.id)]), mod; args...)
#TODO: for more expr, pack them into a tuple

gather(sim::Simulation, expr, mod::Module = AstroNbodySim) = gather(sim.pids, :(registry[$(sim.id)].$expr), mod)
gather(sim::Simulation, expr1, expr2, mod::Module = AstroNbodySim) = gather(sim.pids, :(registry[$(sim.id)].$expr1.$expr2), mod)
gather(sim::Simulation, f::Function, expr, mod::Module = AstroNbodySim; args...) = gather(f, sim.pids, :(registry[$(sim.id)].$expr), mod; args...)
gather(sim::Simulation, f::Function, mod::Module = AstroNbodySim; args...) = gather(f, sim.pids, :(registry[$(sim.id)]), mod; args...)

sum(sim::Simulation, expr, mod::Module = AstroNbodySim) = sum(gather(sim, expr, mod))
sum(sim::Simulation, expr1, expr2, mod::Module = AstroNbodySim) = sum(gather(sim, expr1, expr2, mod))

minimum(sim::Simulation, expr, mod::Module = AstroNbodySim) = minimum(gather(sim, expr, mod))
minimum(sim::Simulation, expr1, expr2, mod::Module = AstroNbodySim) = minimum(gather(sim, expr1, expr2, mod))
maximum(sim::Simulation, expr, mod::Module = AstroNbodySim) = maximum(gather(sim, expr, mod))
maximum(sim::Simulation, expr1, expr2, mod::Module = AstroNbodySim) = maximum(gather(sim, expr1, expr2, mod))

function send_buffer(sim::Simulation)
    # Reduce communication blocking
    # Move myid to last
    src = myid()
    circpids = circshift(sim.pids, length(sim.pids) - findfirst(x->x==src, sim.pids))

    BufferType = typeof(first(sim.buffer.sendbuffer).second)
    for target in circpids[1:end-1]
        sim.recvbuffer[target] = (Distributed.remotecall_eval(AstroNbodySim, target, :(registry[$(sim.id)].buffer.sendbuffer[$src])))::BufferType
    end
end

function send_buffer(sim::Simulation, sendsymbol::Symbol, recvsymbol::Symbol)
    src = myid()
    circpids = circshift(sim.pids, length(sim.pids) - findfirst(x->x==src, sim.pids))

    BufferType = typeof(first(getfield(sim.buffer, sendsymbol)).second)
    for target in circpids[1:end-1]
        remotedata = (Distributed.remotecall_eval(AstroNbodySim, target, :(registry[$(sim.id)].buffer.$(sendsymbol)[$src])))::BufferType
        setindex!(getfield(sim.buffer, recvsymbol), remotedata, target)
    end
end

function send_force_buffer(sim::Simulation)
    send_buffer(sim, :gravToEval, :gravToEvalRecv)
end

function send_force_result_buffer(sim::Simulation)
    send_buffer(sim, :gravResult, :gravResultRecv)
end

function send_potential_buffer(sim::Simulation)
    send_buffer(sim, :potToEval, :potToEvalRecv)
end

function send_potential_result_buffer(sim::Simulation)
    send_buffer(sim, :potResult, :potResultRecv)
end

function clear_local()
    empty!(registry)
end

"""
    function clear(pids = procs())

Clear distributed memories in `AstroNbodySim.registry`, `AstroNbodySim.PhysicalTrees.registry`
"""
function clear(pids = procs())
    bcast(pids, clear_local)
    PhysicalTrees.clear(pids)
end

# Have to overload PhysicalTrees.send_buffer, since the module structure might be Main.AstroNbodySim.PhysicalTrees.
# Remote call might not find Main.PhysicalTrees
#function send_buffer(tree::AbstractTree)
#    # Reduce communication blocking
#    # Move myid to last
#    src = myid()
#    circpids = circshift(tree.pids, length(tree.pids) - findfirst(x->x==src, tree.pids))
#
#    for target in circpids[1:end-1]
#        tree.recvbuffer[target] = Distributed.remotecall_eval(AstroNbodySim.PhysicalTrees, target, :(registry[$(tree.id)].sendbuffer[$src]))
#    end
#end

# So is the same with PhysicalMeshes