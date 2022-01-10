function saverestart(sim::Simulation, ::Union{DirectSum, Tree}, ::CPU)
    if !(sim.id.first in sim.pids)
        @everywhere sim.id.first Main.AstroNbodySim.FileIO.save(
            joinpath($(sim.config.output.dir), "restart.$(Main.AstroNbodySim.Distributed.myid()).jld2"),
            Dict("sim" => Main.AstroNbodySim.registry[$(sim.id)])
        )
    end
    @everywhere sim.pids Main.AstroNbodySim.FileIO.save(
        joinpath($(sim.config.output.dir), "restart.$(Main.AstroNbodySim.Distributed.myid()).jld2"),
        Dict("sim" => Main.AstroNbodySim.registry[$(sim.id)])
    )
end

function saverestart(sim::Simulation, ::DirectSum, ::GPU)
    Main.AstroNbodySim.FileIO.save(
        joinpath(sim.config.output.dir, "restart.$(myid()).jld2"),
        Dict("sim" => sim)
    )
end

function saverestart(sim::Simulation, ::Union{FFT, FDM}, ::CPU)
    Main.AstroNbodySim.FileIO.save(
        joinpath(sim.config.output.dir, "restart.$(myid()).jld2"),
        Dict("sim" => sim)
    )
end


function saverestart(sim::Simulation)
    if sim.visinfo.Realtime
        sim.visinfo.PlotData = nothing
        sim.visinfo.scene = nothing
    end

    saverestart(sim, sim.config.solver.grav, sim.config.device.type)
end


function redistribute_data(sim::Simulation, ::Gravity)
end

function redistribute_data(sim::Simulation, ::Tree)
    # We simply rebuild the tree to avoid errors, and warm up the code
    #! Assume that the holder of sim is also the holder of tree, do not handle errors for other cases
    newtreeid = Main.PhysicalTrees.next_treeid()
    sim.simdata.tree.id = newtreeid

    @everywhere sim.pids PhysicalTrees.registry[$(newtreeid)] = AstroNbodySim.registry[$(sim.id)].tree

    # Check the holder
    if !haskey(Main.PhysicalTrees.registry, newtreeid)
        Core.eval(Main.PhysicalTrees, :(registry[newtreeid] = $(sim.simdata.tree)))
    end

    sim.simdata.tree.pids = sim.pids
    bcast(sim.simdata.tree, :id, :newtreeid)
    bcast(sim.simdata.tree, :pids, sim.pids)
    rebuild(sim.simdata.tree)
end

function loadrestartremote(sim::Simulation, dir::String, pids::Array{Int64,1})
    newid = next_simulationid()

    toload = filter(x->x!=sim.id.first, sim.pids)
    newpids = filter(x->x!=myid(), pids)
    @sync begin
        asyncmap(1:length(toload)) do i
            @info "loading restart file on process $(newpids[i])"
            @everywhere newpids[i] AstroNbodySim.registry[$(newid)] = FileIO.load($(joinpath(dir, "restart.$(toload[i]).jld2")), "sim")
        end
    end

    registry[newid] = sim
    
    sim.id = newid

    bcast(sim, :pids, pids)
    sim.pids = pids

    bcast(sim, :id, newid)

    redistribute_data(sim, sim.config.solver.grav)
    return sim
end

function check_restartable(sim::Simulation, pids::Array{Int64,1};
                           autoload = false)
    @info "Check restart file"
    #if !(sim.isholder)
    #    error("Need to load the holder restart file first!")
    #end

    # return true if modules are loaded
    status = gather(pids, :(isdefined(Main, :PhysicalTrees) && 
                       isdefined(Main, :FileIO) && 
                       isdefined(Main, :AstroNbodySim)))

    for i in 1:length(pids)
        if !status[i]
            if autoload
                @warn "Loading needed modules on process $(pids[i])"
                @everywhere pids[i] Core.eval(Main, :(using AstroNbodySim, PhysicalParticles, FileIO))
            else
                error("module environment on process $(pids[i]) is not satisfied!\n",
                      "Needed modules under Main: AstroNbodySim, PhysicalTrees, FileIO")
            end
        end
    end

    # myid has to be at the first of pids
    if sim.id.first in sim.pids
        if !(myid() in pids)
            @warn "Process $(myid()) has to be at the first of pids"
            pids[1] = myid()
        else
            index = findfirst(x->x==myid(), pids)
            if index != 1
                @warn "Process $(myid()) has to be at the first of pids"
                # switch
                pids[index] = pids[1]
                pids[1] = myid()
            end
        end
    end
end

function loadrestart(filename::String, pids::Array{Int64,1})
    sim = FileIO.load(filename, "sim")
    check_restartable(sim, pids)
    return loadrestartremote(sim, dirname(abspath(filename)), pids)
end

function loadrestart(filename::String)
    sim = FileIO.load(filename, "sim")
    check_restartable(sim, sim.pids)
    return loadrestartremote(sim, dirname(abspath(filename)), sim.pids)
end