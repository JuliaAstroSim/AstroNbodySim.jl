function preprocessdata(sim::Simulation, ::Gravity, ::CPU)
    bcast(sim, preferunits, :config, :units)
end

function preprocessdata(sim::Simulation, ::Gravity, ::GPU)
    
end

function preprocessdata(sim::Simulation, ::DirectSum, ::GPU)
    preferunits(sim.config.units)
end

"""
    function run(sim::Simulation)

This function does all the work for you:
1. `prepare` sim environment and preprocess data
2. evaluate force, initialize timesteps, output
3. `mainloop`
"""
function run(sim::Simulation)
    timeinfo = sim.timeinfo
    config = sim.config

    TimeIntegrationAlgorithm = config.time.algorithm
    GravSolver = config.solver.grav
    Device = config.device.type

    # Things to do before sim
    preprocessdata(sim, sim.config.solver.grav, sim.config.device.type)

    # Plot
    if sim.visinfo.Realtime
        makie_scatter(sim, GravSolver, Device)
    end

    # Output folder
    mkpathIfNotExist(sim.config.output.dir)
    @info "Snapshot are saved to " * sim.config.output.dir

    # Remove 'stop' file
    stopfile = joinpath(sim.config.output.dir, "stop") 
    if isfile(stopfile)
        rm(stopfile)
    end

    # Logging, timing, profiling, analyzing, progress
    setuploggers(sim)

    # First run, warm up
    compute_force(sim, GravSolver, Device)

    init_timesteps(sim, GravSolver, Device)

    if sim.config.output.pot
        compute_potential(sim, GravSolver, Device)
    end

    ## Leapfrog first half forward kick
    if sim.config.time.algorithm isa Leapfrog
        leapfrog_half_kick(sim, +1, GravSolver, Device)
    end

    # Save initial conditions
    sim.config.output.func(sim, sim.config.output.type)

    # Analysers
    if !isempty(sim.loginfo.analysers)
        analysis = []
    end

    # Main loop
    TimeBetweenRestarts = config.output.SaveRestartFreq * (config.time.End - config.time.Begin)
    NextSaveRestart = TimeBetweenRestarts
    state = true
    while state
        t_TOTAL = begin_timer(sim, TOTAL)
        begin_timer(sim, FORCE)
        begin_timer(sim, DRIFT)
        begin_timer(sim, KICK)
        begin_timer(sim, ANALYSIS)
        begin_timer(sim, OUTPUT)
        begin_timer(sim, PLOT)
        
        # force computation, time integration, logging of various types of simulation
        step(sim, GravSolver, Device) 

        # Do whatever you want with the entire data of simulation !!!
        t_ANALYSIS = time_ns()
        if !isempty(sim.loginfo.analysers)
            push!(analysis, write_analysis(sim))
        end
        add_timer(sim, ANALYSIS, t_ANALYSIS, time_ns())
        
        # Plot
        if sim.visinfo.Realtime && time() - sim.visinfo.last_plot_time > sim.visinfo.RenderTime
            t_PLOT = time_ns()
            update_makie_plot(sim, GravSolver, Device)
            add_timer(sim, PLOT, t_PLOT, time_ns())
        end

        # save restart
        if sim.config.output.SaveRestart && sim.timeinfo.system_time_float >= NextSaveRestart
            saverestart(sim)
            NextSaveRestart += TimeBetweenRestarts
        end

        # update progress bar
        if sim.config.loggingmode == ProgressMode()
            update_progress(sim, config.time.step; showvalues = [
                ("Step", timeinfo.stepcount),
                ("System time", timeinfo.system_time_float),
                ("Timestep", timeinfo.system_time_float - timeinfo.last_system_time_float)
            ])
        end

        # interrupt
        if need_to_interrupt(sim.config.output.dir)
            @warn "\n\n\n\n\nInterrupted by user"
            state = false
            closeloggers(sim)
        end
        if timeinfo.system_time_float >= config.time.End
            state = false
        end

        # Update system time and counter
        timeinfo.stepcount += 1
        timeinfo.last_system_time_float = timeinfo.system_time_float

        add_timer(sim, TOTAL, t_TOTAL, time_ns())
        write_timing(sim)
    end

    ## Leapfrog last half backward kick
    if sim.config.time.algorithm isa Leapfrog
        leapfrog_half_kick(sim, -1, GravSolver, Device)
    end

    # Finish simulation
    closeloggers(sim)
    println()
    println()
    #@info "Simulation finished"

    if !isempty(sim.loginfo.analysers)
        return collect(hcat(analysis...)')
    end
end

function restart(sim::Simulation)
    @info "Restarting run"
    restartloggers(sim)
    mainloop(sim)
end