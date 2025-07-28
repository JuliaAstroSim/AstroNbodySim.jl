"""
$(TYPEDSIGNATURES)

Compute accelerations by second-order central differencing of potential on the mesh.
Results are stored in `acc` and return `nothing`.
"""
function compute_acc(m::MeshCartesianStatic, dim::Val{3}, mode::VertexMode)
    config = m.config

    # acc on the boundaries are omitted. This does no harm because we have enlarged the mesh extent.
    m.acc.x .= -diff_central_x(config.Δ[1], m.phi)
    m.acc.y .= -diff_central_y(config.Δ[2], m.phi)
    m.acc.z .= -diff_central_z(config.Δ[3], m.phi)
    return nothing
end

function compute_acc(m::MeshCartesianStatic)
    compute_acc(m, Val(m.config.dim), m.config.mode)
end

function compute_force(sim::Simulation, GravSolver::Union{FDM, FFT, ML}, Device::DeviceType)
    t_FORCE = time_ns()

    G = sim.config.constants.G
    ACC0 = sim.config.constants.ACC0
    m = sim.simdata
    data = m.data

    if GravSolver isa FDM
        fdm_poisson(m, Val(m.config.dim), G, m.config.mode, Device, m.config.boundary, sim.config.grav.sparse)
    elseif GravSolver isa FFT
        fft_poisson(m, G)
    elseif GravSolver isa ML
        cnn_poisson(m, sim.config.solver.data.u, sim.config.solver.data)
    end
    
    compute_acc(m)

    #TODO background force

    # Solve QUMOND on mesh
    if sim.config.grav.model isa QUMOND
        QUMOND_acc!(m, ACC0, G)
    end

    # Assign acc to inbound particles
    assignparticle(m.data, m, :Acc, :acc)

    # Apply MOND to all particles
    if sim.config.grav.model isa MOND1983Milgrom
        if Device isa GPU
            nuIndex = sim.config.grav.MOND_nuIndex
            NumThreads = sim.config.device.GPU_NumThreads
            NumBlocks = gpu_blocks(sim.config, length(data))
            @cuda threads=NumThreads blocks=NumBlocks mond_Milgrom1983_gpu_kernel(data, nuIndex, ACC0)
        elseif Device isa CPU
            mond_Milgrom1983(sim, data)
        end
    end

    # When use DS to handle outbound particles, if QUMOND, we apply Milgrom 1983 formula to outbound particles
    if sim.config.grav.model isa QUMOND && sim.config.grav.outbound isa DS
        list = outbound_list(m)
        nuIndex = sim.config.grav.MOND_nuIndex
        
        if Device isa GPU
            NumThreads = sim.config.device.GPU_NumThreads
            NumBlocks = gpu_blocks(sim.config, length(data))
            @cuda threads=NumThreads blocks=NumBlocks mond_Milgrom1983_gpu_kernel(data, nuIndex, ACC0)
        elseif Device isa CPU
            Threads.@threads for k in 1:Threads.nthreads()
                Head, Tail = split_block(length(list), k, Threads.nthreads())
                for i in Head:Tail
                    @inbounds NewtonAcc = data.Acc[list[i]]
                    @inbounds data.Acc[list[i]] = NewtonAcc * nu(norm(NewtonAcc) / ACC0, nuIndex)
                end
            end
        end
    end

    add_timer(sim, "FORCE", t_FORCE, time_ns())
end

function compute_force(sim::Simulation, pos::Array{T,N}, SoftLength::Number, GravSolver::Union{FFT,FDM}, Device::CPU) where T<:AbstractPoint3D where N
    return mesh2particle.(Ref(sim.simdata.pos), sim.simdata.config, Ref(sim.simdata.acc), pos, sim.simdata.config.mode, sim.simdata.config.assignment)
end