function drift_particles_adaptive_timestep(sim::Simulation, next_int::Int , GravSolver::Union{FFT, FDM}, Device::CPU)
    timeinfo = sim.timeinfo
    timeinfo.last_system_time_int = timeinfo.system_time_int
    timeinfo.system_time_int = next_int
    timeinfo.system_time_float = sim.config.time.Begin + timeinfo.system_time_int * sim.config.time.step.TimeInterval
    drift_particles_kernel(sim, GravSolver, Device)
end

function drift_particles_const_timestep(sim::Simulation, next_float::Number , GravSolver::Union{FFT, FDM}, Device::CPU)
    timeinfo = sim.timeinfo
    timeinfo.last_system_time_float = timeinfo.system_time_float
    timeinfo.system_time_float = next_float
    drift_particles_kernel(sim, GravSolver, Device)
end

function outbound_rule(sim::Simulation, m::MeshCartesianStatic, ::Delete)
    list = outbound_list(m)
    StructArrays.foreachfield(v->deleteat!(v, list), sim.simdata.data)
end


function outbound_rule(sim::Simulation, m::MeshCartesianStatic, ::DS)
    list = outbound_list(m)
    data = m.data
    G = sim.config.constants.G

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(list), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Acc[list[i]] = compute_force_at_point(data.Pos[list[i]], data, G, softlen(data.Collection[list[i]], sim))
        end
    end
end


function outbound_rule(sim::Simulation, m::MeshCartesianStatic, ::CoarseMesh)
    G = sim.config.constants.G
    ACC0 = sim.config.constants.ACC0
    list = outbound_list(m)
    device = sim.config.device.type

    #TODO check same setup with m
    coarsemesh = MeshCartesianStatic(m.data;
        mode = m.config.mode,
        assignment = m.config.assignment,
        boundary = m.config.boundary,
        Nx = 10,
        Ny = 10,
        Nz = 10,
        assign = true,
        gpu = device isa GPU ? true : false,
        enlarge = 1.2,
    )

    if GravSolver isa FDM
        fdm_poisson(coarsemesh, Val(coarsemesh.config.dim), G, coarsemesh.config.mode, device, coarsemesh.config.boundary, sim.config.grav.sparse)
    elseif GravSolver isa FFT
        fft_poisson(coarsemesh, G)
    end

    compute_acc(coarsemesh, Val(coarsemesh.config.dim), coarsemesh.config.mode)

    #TODO background force

    if sim.config.grav.model isa QUMOND
        QUMOND_acc!(coarsemesh, ACC0, G)
    end

    # Assign outbound particles
    #TODO this is unsafe assign
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(list), k, Threads.nthreads())
        for i in Head:Tail
            pos = SVector(m.data.Pos[list[i]])
            getproperty(m.data, :Acc)[list[i]] = mesh2particle(coarsemesh, pos, :acc, coarsemesh.config.mode, coarsemesh.config.assignment)
        end
    end

    #TODO GPU
end

function step(sim::Simulation, GravSolver::Union{FFT, FDM}, Device::DeviceType)
    m = sim.simdata

    # Drift and output
    find_next_sync_point_and_drift(sim, sim.config.time.step, GravSolver, Device)

    # Outbound
    outbound_rule(sim, m, sim.config.grav.outbound)

    # Reassign mesh
    #TODO move mesh.rho
    assignmesh(m)

    compute_force(sim, GravSolver, Device)

    # Kick
    advance_and_find_timestep(sim, GravSolver, Device)
end