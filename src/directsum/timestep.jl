function drift_particles_adaptive_timestep(sim::Simulation, next_int::Int64, GravSolver::DirectSum, Device::CPU)
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

    # Plot
    if sim.visinfo.Realtime && time() - sim.visinfo.last_plot_time > sim.visinfo.RenderTime
        update_makie_plot(sim, GravSolver, Device)
    end
end

function drift_particles_const_timestep(sim::Simulation, next_float::Number, GravSolver::DirectSum, Device::CPU)
    timeinfo = sim.timeinfo

    # Update system time
    timeinfo.last_system_time_float = timeinfo.system_time_float
    bcast(sim, :timeinfo, :last_system_time_float, timeinfo.system_time_float)

    timeinfo.system_time_float = next_float
    bcast(sim, :timeinfo, :system_time_float, next_float)

    # Drift
    bcast(sim, drift_particles_kernel, args = (GravSolver, Device))

    # Plot
    if sim.visinfo.Realtime && time() - sim.visinfo.last_plot_time > sim.visinfo.RenderTime
        update_makie_plot(sim, GravSolver, Device)
    end
end


function step(sim::Simulation, GravSolver::DirectSum, Device::CPU)
    #timeinfo = sim.timeinfo
    #write_log(sim, string(
    #    "\nStep ", timeinfo.stepcount, "    System time: ", timeinfo.system_time_float
    #))

    # Drift and output
    find_next_sync_point_and_drift(sim, sim.config.time.step, GravSolver, Device)

    compute_force(sim, GravSolver, Device)

    # Kick
    advance_and_find_timestep(sim, GravSolver, Device)
end