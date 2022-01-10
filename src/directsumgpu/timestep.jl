function get_timestep_gpu(acc::AbstractPoint, SofteningLength::Number, TimestepConfig::AdaptiveTimestep)
    a2 = norm(acc)
    dt = sqrt(2 * TimestepConfig.ErrTolTimestep * SofteningLength / a2)

    if dt > TimestepConfig.MaxStep
        dt = TimestepConfig.MaxStep
    end

    #! cannot print quantity
    if dt < TimestepConfig.MinStep
        @cuprintln("Timestep < MinStep")
        #@cuprintln("Timestep < MinStep = ", TimestepConfig.MinStep, ", particle = ", p, ", dt = ", dt)
    end

    ti_step = dt / TimestepConfig.TimeInterval
    if ti_step <= 0 || ti_step > TimestepConfig.TimeBase
        @cuprintln("Wrong integer timestep")
        #@cuprintln("Wrong integer timestep: ", ti_step, ", p = ", p)
    end

    ti_min = TimestepConfig.TimeBase
    while ti_min > ti_step
        ti_min >>= 1
    end
    return ti_min
end

function get_timestep_gpu(acc::AbstractPoint, SofteningLength::Number, TimestepConfig::ConstantTimestep)
    return 0
end

function init_timesteps_gpu_kernel(particles::StructArray, ForceSofteningTable::CuDeviceArray, TimestepConfig::AbstractTimestepConfig)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Ti_endstep[i] = get_timestep_gpu(particles.Acc[i], ForceSofteningTable[Int(particles.Collection[i])], TimestepConfig)
    end
    return nothing
end

function init_timesteps(sim::Simulation, ::Gravity, ::GPU)
    dataGPU = get_local_data(sim)
    NumThreads = sim.config.device.GPU_NumThreads
    cuForceTable = cu(sim.config.grav.ForceSofteningTable)

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks init_timesteps_gpu_kernel(dataGPU, cuForceTable, sim.config.time.step)
end

function drift_particles_gpu_kernel(particles::StructArray, dt::Number)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Pos[i] += particles.Vel[i] * dt
    end
    return nothing
end

function drift_particles_kernel(sim::Simulation, TimestepConfig::AdaptiveTimestep, ::Gravity, ::GPU)
    timeinfo = sim.timeinfo
    dataGPU = get_local_data(sim)
    dt = timeinfo.dt = (timeinfo.system_time_int - timeinfo.last_system_time_int) * TimestepConfig.TimeInterval
    NumThreads = sim.config.device.GPU_NumThreads

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks drift_particles_gpu_kernel(dataGPU, dt)
end

function drift_particles_kernel(sim::Simulation, TimestepConfig::ConstantTimestep, ::Gravity, ::GPU)
    dt = TimestepConfig.dt
    dataGPU = get_local_data(sim)
    NumThreads = sim.config.device.GPU_NumThreads

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks drift_particles_gpu_kernel(dataGPU, dt)
end

function drift_particles_adaptive_timestep(sim::Simulation, next_int::Integer, GravSolver::Gravity, Device::GPU)
    timeinfo = sim.timeinfo
    timeinfo.last_system_time_int = timeinfo.system_time_int
    timeinfo.system_time_int = next_int
    timeinfo.system_time_float = sim.config.time.Begin + timeinfo.system_time_int * sim.config.time.step.TimeInterval

    drift_particles_kernel(sim, sim.config.time.step, GravSolver, Device)
end

function drift_particles_const_timestep(sim::Simulation, next_float::Number, GravSolver::Gravity, Device::GPU)
    timeinfo = sim.timeinfo
    timeinfo.last_system_time_float = timeinfo.system_time_float
    timeinfo.system_time_float = next_float

    drift_particles_kernel(sim, sim.config.time.step, GravSolver, Device)
end

#function copy_endstep_to_gpu(particles::CuDeviceArray, endsteps::CuDeviceArray)
#    i = threadIdx().x
#    endsteps[i] = particles[i].Ti_endstep
#    return nothing
#end

function get_endstep(p::AbstractParticle)
    return p.Ti_endstep
end

function find_min_endstep(sim::Simulation, ::Gravity, ::GPU)
    min_endstep = sim.config.time.step.TimeBase
    dataGPU = get_local_data(sim)
    endsteps = dataGPU.Ti_endstep
    m = minimum(endsteps)
    if min_endstep > m
        min_endstep = m
    end

    sim.timeinfo.min_endstep = min_endstep
    return min_endstep
end

function check_sync_point(sim::Simulation, ::Gravity, ::GPU)
    check_sync_point_local(sim)
end



function advance_and_find_timestep_gpu_kernel(particles::StructArray, ForceSofteningTable::CuDeviceArray, TimestepConfig::AdaptiveTimestep, system_time_int::Integer)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles) && particles.Ti_endstep[i] == system_time_int
        ti_step = get_timestep_gpu(particles.Acc[i], ForceSofteningTable[Int(particles.Collection[i])], TimestepConfig)

        if system_time_int == TimestepConfig.TimeBase
            ti_step = 0
        end

        if TimestepConfig.TimeBase - system_time_int < ti_step
            ti_step = TimestepConfig.TimeBase - system_time_int
        end

        tstart = (particles.Ti_begstep[i] + particles.Ti_endstep[i]) / 2
        tend = particles.Ti_endstep[i] + 0.5 * ti_step

        dt = (tend - tstart) * TimestepConfig.TimeInterval
        
        # Kick
        particles.Vel[i] += particles.Acc[i] * dt
        particles.Ti_begstep[i] = particles.Ti_endstep[i] #! do not update these two in the same time
        particles.Ti_endstep[i] += ti_step
    end
    return nothing
end

function advance_and_find_timestep(sim::Simulation, dataGPU::StructArray, TimestepConfig::AdaptiveTimestep, ::Gravity, ::GPU)
    NumThreads = sim.config.device.GPU_NumThreads
    cuForceTable = cu(sim.config.grav.ForceSofteningTable)

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks advance_and_find_timestep_gpu_kernel(dataGPU, cuForceTable, TimestepConfig, sim.timeinfo.system_time_int)
end

function advance_and_find_timestep_gpu_kernel(particles::StructArray, dt::Number)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Vel[i] += particles.Acc[i] * dt
    end
    return nothing
end

function advance_and_find_timestep(sim::Simulation, dataGPU::StructArray, TimestepConfig::ConstantTimestep, ::Gravity, ::GPU)
    dt = sim.timeinfo.dt = TimestepConfig.dt
    NumThreads = sim.config.device.GPU_NumThreads

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks advance_and_find_timestep_gpu_kernel(dataGPU, dt)
end

function advance_and_find_timestep(sim::Simulation, GravSolver::Gravity, Device::GPU)
    advance_and_find_timestep(sim, get_local_data(sim), sim.config.time.step, GravSolver, Device)
end



function leapfrog_half_kick_gpu_kernel(particles::StructArray, ForceSofteningTable::CuDeviceArray, TimestepConfig::AdaptiveTimestep, system_time_int::Integer, forward::Int)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Vel[i] += particles.Acc[i] * 0.5 * forward * (particles.Ti_endstep[i] - particles.Ti_begstep[i]) * TimestepConfig.TimeInterval
    end
    return nothing
end

function leapfrog_half_kick_local(sim::Simulation, dataGPU::StructArray, TimestepConfig::AdaptiveTimestep, forward::Int, ::Gravity, ::GPU)
    NumThreads = sim.config.device.GPU_NumThreads
    cuForceTable = cu(sim.config.grav.ForceSofteningTable)

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks leapfrog_half_kick_gpu_kernel(dataGPU, cuForceTable, TimestepConfig, sim.timeinfo.system_time_int, forward)
end

function leapfrog_half_kick_gpu_kernel(particles::StructArray, dt::Number, forward::Int)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Vel[i] += particles.Acc[i] * 0.5 * forward * dt
    end
    return nothing
end

function leapfrog_half_kick_local(sim::Simulation, dataGPU::StructArray, TimestepConfig::ConstantTimestep, forward::Int, ::Gravity, ::GPU)
    dt = sim.timeinfo.dt = TimestepConfig.dt
    NumThreads = sim.config.device.GPU_NumThreads

    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    @cuda threads=NumThreads blocks=NumBlocks leapfrog_half_kick_gpu_kernel(dataGPU, dt, forward)
end

function leapfrog_half_kick(sim::Simulation, forward::Int, GravSolver::Gravity, Device::GPU)
    leapfrog_half_kick_local(sim, get_local_data(sim), sim.config.time.step, forward, GravSolver, Device)
end



function step(sim::Simulation, GravSolver::DirectSum, Device::GPU)
    #TODO logging
    # Drift and output
    find_next_sync_point_and_drift(sim, sim.config.time.step, GravSolver, Device)

    compute_force(sim, GravSolver, Device)

    # Kick
    advance_and_find_timestep(sim, GravSolver, Device)
end