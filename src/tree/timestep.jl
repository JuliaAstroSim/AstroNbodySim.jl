#function need_to_rebuild_tree(sim::Simulation)
#    NumForceUpdateTotal = sum(sim, :physics, :NumForceUpdateSinceLast)
#    if NumForceUpdateTotal > sim.simdata.tree.mutable.NumTotal * sim.DomainUpdateFrequency
#        return true
#    else
#        return false
#    end
#end

function rebuild_tree(sim::Simulation)
    bcast(sim, :physics, :NumForceUpdateSinceLast, 0)
    rebuild(sim.simdata.tree)
end

function drift_particles_adaptive_timestep(sim::Simulation, next_int::Int64, GravSolver::Tree, Device::CPU)
    timeinfo = sim.timeinfo

    # Update system time
    timeinfo.last_system_time_int = timeinfo.system_time_int
    bcast(sim, :timeinfo, :last_system_time_int, timeinfo.system_time_int)

    timeinfo.system_time_int = next_int
    bcast(sim, :timeinfo, :system_time_int, next_int)

    timeinfo.system_time_float = sim.config.time.Begin + timeinfo.system_time_int * sim.config.time.step.TimeInterval
    bcast(sim, :timeinfo, :system_time_float, timeinfo.system_time_float)

    # Drift
    bcast(sim, drift_particles_kernel, args = (GravSolver, Device))

    rebuild_tree(sim)
end

function drift_particles_const_timestep(sim::Simulation, next_float::Number, GravSolver::Tree, Device::CPU)
    timeinfo = sim.timeinfo

    # Update system time
    timeinfo.last_system_time_float = timeinfo.system_time_float
    bcast(sim, :timeinfo, :last_system_time_float, timeinfo.system_time_float)

    timeinfo.system_time_float = next_float
    bcast(sim, :timeinfo, :system_time_float, next_float)

    # Drift
    bcast(sim, drift_particles_kernel, args = (GravSolver, Device))

    rebuild_tree(sim)
end



function step(sim::Simulation, GravSolver::Tree, Device::CPU)
    #TODO logging
    # Drift and output
    find_next_sync_point_and_drift(sim, sim.config.time.step, GravSolver, Device)

    compute_force(sim, GravSolver, Device)

    # Kick
    advance_and_find_timestep(sim, GravSolver, Device)
end