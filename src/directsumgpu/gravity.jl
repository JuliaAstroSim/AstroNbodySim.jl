function gpu_threads(config::SimConfig)
    return config.device.GPU_NumThreads
end

function gpu_threads(sim::Simulation)
    return gpu_threads(sim.config)
end


function gpu_blocks(NumThreads::Int, NumTotal::Int)
    return div(NumTotal, NumThreads) + 1
end

function gpu_blocks(config::SimConfig, NumTotal::Int)
    return div(NumTotal, gpu_threads(config)) + 1
end

function gpu_blocks(sim::Simulation)
    return gpu_blocks(sim.config, length(sim.simdata))
end

function compute_force(sim::Simulation, Algorithm::Tiled)
    
end

Base.@propagate_inbounds function compute_force_gpu_kernel(particles::StructVector, sources::StructVector, G::Number, ForceSofteningTable::CuDeviceVector)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        a = particles.Acc[i]
        pos = particles.Pos[i]
        SoftLen = ForceSofteningTable[Int(particles.Collection[i])]
        for k in eachindex(sources)
            a += compute_force_from_particle(pos, sources.Pos[k], sources.Mass[k], G, SoftLen)
        end
        particles.Acc[i] = a
    end
    return nothing
end

function compute_force(sim::Simulation, Algorithm::AllPairs)
    dataGPU = sim.simdata
    G = sim.config.constants.G
    cuForceTable = cu(sim.config.grav.ForceSofteningTable)

    NumThreads = sim.config.device.GPU_NumThreads
    NumBlocks = gpu_blocks(sim.config, length(dataGPU))
    dataGPU.Acc .*= 0.0
    @cuda threads=NumThreads blocks=NumBlocks compute_force_gpu_kernel(dataGPU, dataGPU, G, cuForceTable)
end

function apply_background_force_gpu_kernel(f::Function, particles::StructVector)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        particles.Acc[i] += f(particles.Pos[i])
    end
    return nothing
end

function mond_Milgrom1983_gpu_kernel(particles::StructVector, nuIndex, ACC0)
    i = threadIdx().x + blockDim().x * (blockIdx().x - 1)
    if i <= length(particles)
        NewtonAcc = particles[i].Acc
        particles.Acc[i] = NewtonAcc * nu(norm(NewtonAcc) / ACC0, nuIndex)
    end
    return nothing
end

function postprocessing_force(sim::Simulation, GravSolver::DirectSum, Device::GPU)
    NumThreads = sim.config.device.GPU_NumThreads
    NumBlocks = gpu_blocks(sim.config, length(sim.simdata))

    #TODO: unsupported dynamic function invocation
    #! bgforce function should not use global variables
    if sim.config.grav.model isa MOND1983Milgrom
        nuIndex = sim.config.grav.MOND_nuIndex
        ACC0 = sim.config.constants.ACC0
        @cuda threads=NumThreads blocks=NumBlocks mond_Milgrom1983_gpu_kernel(sim.simdata, nuIndex, ACC0)
    end
end

function compute_force(sim::Simulation, GravSolver::DirectSum, Device::GPU)
    t_FORCE = time_ns()
    compute_force(sim, sim.config.device.GPU_Algorithm)
    postprocessing_force(sim, GravSolver, Device)
    add_timer(sim, FORCE, t_FORCE, time_ns())
end