function get_timestep(acc::AbstractPoint, SofteningLength::Number, TimestepConfig::AdaptiveTimestep)
    ac = norm(acc)
    dt = sqrt(2 * TimestepConfig.ErrTolTimestep * SofteningLength / ac) #TODO: what if SofteningLength == 0?

    if dt > TimestepConfig.MaxStep
        dt = TimestepConfig.MaxStep
    end

    if dt < TimestepConfig.MinStep
        error("Timestep < MinStep = ", TimestepConfig.MinStep, #=", particle = ", p,=# ", dt = ", dt)
    end

    ti_step = dt / TimestepConfig.TimeInterval
    if ti_step <= 0 #|| ti_step > TimestepConfig.TimeBase
        error("Wrong integer timestep: ", ti_step, #=", p = ", p=#)
    end

    ti_min = TimestepConfig.TimeBase
    while ti_min > ti_step
        ti_min >>= 1
    end
    return ti_min
end

function get_timestep(acc::AbstractPoint, SofteningLength::Number, TimestepConfig::ConstantTimestep)
    return 0
end

function init_timesteps_local(sim::Simulation)
    data = get_local_data(sim)

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Ti_endstep[i] = get_timestep(data.Acc[i], softlen(data.Collection[i], sim), sim.config.time.step)
        end
    end
    sim.physics.NumForceUpdateSinceLast = 0
end

function init_timesteps(sim::Simulation, ::Gravity, ::CPU)
    bcast(sim, init_timesteps_local)
end
init_timesteps(sim::Simulation) = init_timesteps(sim, sim.config.solver.grav, sim.config.device.type)


Base.@propagate_inbounds function advance_and_find_timestep_local(sim::Simulation, data::StructArray, TimestepConfig::AdaptiveTimestep, ::Gravity, ::CPU)
    TimeBase = TimestepConfig.TimeBase
    TimeInterval = TimestepConfig.TimeInterval
    timeinfo = sim.timeinfo

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            if data.Ti_endstep[i] == timeinfo.system_time_int
                ti_step = get_timestep(data.Acc[i], softlen(data.Collection[i], sim), TimestepConfig)

                if timeinfo.system_time_int == TimeBase
                    ti_step = 0
                end

                if TimeBase - timeinfo.system_time_int < ti_step
                    ti_step = TimeBase - timeinfo.system_time_int
                end

                tstart = (data.Ti_begstep[i] + data.Ti_endstep[i]) / 2
                tend = data.Ti_endstep[i] + 0.5 * ti_step

                dt = (tend - tstart) * TimeInterval
                
                # Kick
                Vel = data.Vel[i] + data.Acc[i] * dt
                Ti_begstep = data.Ti_endstep[i] #! do not update these two in the same time
                Ti_endstep = Ti_begstep + ti_step

                data.Vel[i] = Vel
                data.Ti_begstep[i] = Ti_begstep
                data.Ti_endstep[i] = Ti_endstep
            end
        end
    end
end

function advance_and_find_timestep_local(sim::Simulation, data::StructArray, TimestepConfig::ConstantTimestep, ::Gravity, ::CPU)
    dt = sim.timeinfo.dt = TimestepConfig.dt
    
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Vel[i] += data.Acc[i] * dt
        end
    end
end

function advance_and_find_timestep_local(sim::Simulation, GravSolver::Gravity, Device::CPU)
    advance_and_find_timestep_local(sim, get_local_data(sim), sim.config.time.step, GravSolver, Device)
end

function advance_and_find_timestep(sim::Simulation, GravSolver::Gravity, Device::CPU)
    bcast(sim, advance_and_find_timestep_local, args = (GravSolver, Device))
end



function leapfrog_half_kick_local(sim::Simulation, data::StructArray, TimestepConfig::AdaptiveTimestep, forward::Int, ::Gravity, ::CPU)
    TimeInterval = TimestepConfig.TimeInterval

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Vel[i] += data.Acc[i] * forward * 0.5 * (data.Ti_endstep[i] - data.Ti_begstep[i]) * TimeInterval
        end
    end
end

function leapfrog_half_kick_local(sim::Simulation, data::StructArray, TimestepConfig::ConstantTimestep, forward::Int, ::Gravity, ::CPU)
    dt = sim.timeinfo.dt = TimestepConfig.dt
    
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Vel[i] += data.Acc[i] * forward * 0.5 * dt
        end
    end
end

function leapfrog_half_kick_local(sim::Simulation, forward::Int, GravSolver::Gravity, Device::CPU)
    leapfrog_half_kick_local(sim, get_local_data(sim), sim.config.time.step, forward, GravSolver, Device)
end

function leapfrog_half_kick(sim::Simulation, forward::Int, GravSolver::Gravity, Device::CPU)
    bcast(sim, leapfrog_half_kick_local, args = (forward, GravSolver, Device))
end

function find_min_endstep_local(sim::Simulation)
    min = sim.config.time.step.TimeBase
    data = get_local_data(sim)
    for t in data.Ti_endstep
        if min > t
            min = t
        end
    end
    sim.timeinfo.min_endstep = min
end

function check_sync_point_local(sim::Simulation)
    min_endstep = sim.timeinfo.min_endstep
    data = get_local_data(sim)

    flags = data.Ti_endstep .> min_endstep

    if sum(flags) == length(data)
        sim.outputinfo.syncflag = true
    else
        sim.outputinfo.syncflag = false
    end
end

function check_sync_point(sim::Simulation, GravSolver::Gravity, Device::CPU)
    bcast(sim, check_sync_point_local)
    flags = gather(sim, :outputinfo, :syncflag)
    if sum(flags) > 0
        return sim.outputinfo.syncflag = true
    else
        return sim.outputinfo.syncflag = false
    end
end

function drift_particles_kernel(sim::Simulation, data::StructArray, TimestepConfig::AdaptiveTimestep, ::Gravity, ::CPU)
    timeinfo = sim.timeinfo
    dt = timeinfo.dt = (timeinfo.system_time_int - timeinfo.last_system_time_int) * TimestepConfig.TimeInterval
    
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Pos[i] += data.Vel[i] * dt
        end
    end
end

function drift_particles_kernel(sim::Simulation, data::StructArray, TimestepConfig::ConstantTimestep, ::Gravity, ::CPU)
    #TODO dt should be computed from system_time_float
    dt = sim.timeinfo.system_time_float - sim.timeinfo.last_system_time_float
    
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Pos[i] += data.Vel[i] * dt
        end
    end
end

function drift_particles_kernel(sim::Simulation, GravSolver::Gravity, Device::CPU)
    drift_particles_kernel(sim, get_local_data(sim), sim.config.time.step, GravSolver, Device)
end

function find_next_sync_point_and_drift(sim::Simulation, TimestepConfig::ConstantTimestep, GravSolver::Gravity, Device::DeviceType)
    next_system_time = sim.timeinfo.system_time_float + TimestepConfig.dt

    while sim.timeinfo.next_output_time_float <= next_system_time
        t_DRIFT = time_ns()
        drift_particles_const_timestep(sim, sim.timeinfo.next_output_time_float, GravSolver, Device)
        add_timer(sim, DRIFT, t_DRIFT, time_ns())

        t_OUTPUT = time_ns()
        sim.config.output.func(sim, sim.config.output.type)
        add_timer(sim, OUTPUT, t_OUTPUT, time_ns())

        t_ANALYSIS = time_ns()
        write_analysis(sim)
        add_timer(sim, ANALYSIS, t_ANALYSIS, time_ns())

        sim.timeinfo.next_output_time_float += sim.config.time.BetweenSnapshots
    end

    t_DRIFT = time_ns()
    drift_particles_const_timestep(sim, next_system_time, GravSolver, Device)
    add_timer(sim, DRIFT, t_DRIFT, time_ns())
end

function find_min_endstep(sim::Simulation, GravSolver::Gravity, Device::CPU)
    bcast(sim, find_min_endstep_local)
    min_endstep = minimum(sim, :timeinfo, :min_endstep)
    bcast(sim, :timeinfo, :min_endstep, min_endstep)
    return min_endstep
end

function find_next_sync_point_and_drift(sim::Simulation, TimestepConfig::AdaptiveTimestep, GravSolver::Gravity, Device::DeviceType)
    t_DRIFT = time_ns()
    min_endstep = find_min_endstep(sim, GravSolver, Device)
    check_sync_point(sim, GravSolver, Device)
    add_timer(sim, DRIFT, t_DRIFT, time_ns())

    while sim.timeinfo.next_output_time_int <= min_endstep
        t_DRIFT = time_ns()
        drift_particles_adaptive_timestep(sim, sim.timeinfo.next_output_time_int, GravSolver, Device)
        add_timer(sim, DRIFT, t_DRIFT, time_ns())

        t_OUTPUT = time_ns()
        sim.config.output.func(sim, sim.config.output.type)
        add_timer(sim, OUTPUT, t_OUTPUT, time_ns())

        t_ANALYSIS = time_ns()
        write_analysis(sim)
        add_timer(sim, ANALYSIS, t_ANALYSIS, time_ns())

        sim.timeinfo.next_output_time_int += sim.config.time.BetweenSnapshotsInt
    end

    t_DRIFT = time_ns()
    drift_particles_adaptive_timestep(sim, min_endstep, GravSolver, Device)
    add_timer(sim, DRIFT, t_DRIFT, time_ns())
end