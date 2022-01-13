"""
    tabstring(s::String, tab::String = "    ")

Retract all lines in string `s` by `tab`

"""
function tabstring(s::String, tab::String = "    ")
    return tab * replace(s, "\n" => "\n" * tab)
end

##### OutputConfig #####

"""
$(TYPEDEF)

$(TYPEDFIELDS)

## Notes on Keywords

- Keywords are compatible with `SimConfig`, and can also be overrided by field names of OutputConfig

## Example
```julia
OutputConfig(; dir = "test", type = jld2(), format2 = false, acc = true, pot = true, SaveRestart = true, SaveRestartFreq = 0.33)
```
"""
struct OutputConfig{OT}
    """Output directory. Default: `joinpath(pwd(), "output")`"""
    dir::String
    "Output function. Default: `AstroNbodySim.output`"
    func::Function
    "Output Type: traits of `gadget2`, `hdf5`, `jld2` defined in `AstroIO`. Default: `gadget2()`"
    type::OT

    "If `true`, write in format2 mode for Gadget2 format. Default: `true`"
    format2::Bool
    "If `true`, output acceleration in snapshots. Default: `false`"
    acc::Bool
    "If `true`, output potential in snapshots. Default: `false`"
    pot::Bool

    "If `true`, save restart files regularly. Default: `true`"
    SaveRestart::Bool
    "Regularly after every fraction of the whole simulation time, save restart files. Default: `1.0`"
    SaveRestartFreq::Float64
end

"""
$(TYPEDSIGNATURES)
"""
function OutputConfig(;
    # Keywords from SimConfig
    OutputDir::AbstractString = joinpath(pwd(), "output"),
    OutputFunction::Function = output,
    OutputType::AbstractOutputType = gadget2(),
    format2::Bool = true,
    acc::Bool = false,
    pot::Bool = false,
    SaveRestart::Bool = false,
    SaveRestartFreq::Float64 = 1.0,

    # Keywords to override
    dir::AbstractString = OutputDir,
    func::Function = OutputFunction,
    type::AbstractOutputType = OutputType,
)
    return OutputConfig(
        dir, func, type,
        format2, acc, pot,
        SaveRestart, SaveRestartFreq,
    )
end

function Base.show(io::IO, config::OutputConfig)
    print(io,
        """
        Output config:
                   Directory: $(config.dir)
                    Function: $(config.func)
                        Type: $(traitstring(config.type))
                     Format2: $(config.format2)
                  Output Acc: $(config.acc)
                  Output Pot: $(config.pot)
          Save Restart files: $(config.SaveRestart)
           Save Restart Freq: $(config.SaveRestartFreq)
        """
    )
end



##### GravityConfig #####

"""
$(TYPEDEF)

$(TYPEDFIELDS)

# Example
```julia
GravityConfig(uAstro)
GravityConfig(nothing)
GravityConfig(;ForceSofteningTable=[0.01u"kpc" for i in 1:6])
```
"""
struct GravityConfig{Len, GM, OBL, BC}
    "Gravitational softening lengths for different type of particles"
    ForceSofteningTable::MVector{6,Len}

    "Gravity model. Supported: `Newton`, `MOND1983Milgrom`, `QUMOND`"
    model::GM
    
    "Index of Mondian interpolation function"
    MOND_nuIndex::Float64

    # Mesh
    """
    `::OutboundLimiter`. Choose how to handle particles out of the non-periodic simulation box.
    Supported:
        - `Delete`: delete outbound particles
        - `DS`: use direct summation method to compute forces 
        - `CoarseMesh`: Construct a coarse mesh to overlap all particles
    Default is `DS`.
    """
    outbound::OBL
    "Enlarge the simulation box by `enlarge` compared to the extent of system. Default is `3.0`"
    enlarge::Float64
    "Boundary condition. Supported: `Vacuum`, `Periodic`, `Dirichlet`. Default is `Vacuum`"
    boundary::BC
    "Use sparse matrix to improve performance"
    sparse::Bool
end

"""
$(TYPEDSIGNATURES)
"""
function GravityConfig(units = uAstro;
    # Keywords from SimConfig
    ForceSofteningTable = MVector{6}([
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
    ]),
    MOND_nuIndex = 2.0,
    GravityModel::GravityModel = Newton(),
    OutboundLimiter::OutboundLimiter = DS(),
    EnlargeMesh::Float64 = 2.01,
    BoundaryCondition::BoundaryCondition = Vacuum(),
    sparse::Bool = true,

    # Keywords to override
    model = GravityModel,
    outbound = OutboundLimiter,
    enlarge = EnlargeMesh,
    boundary = BoundaryCondition,
)
    return GravityConfig(
        MVector{6}(ForceSofteningTable),
        model,
        MOND_nuIndex,
        outbound,
        enlarge,
        boundary,
        sparse
    )
end

function Base.show(io::IO, config::GravityConfig)
    print(io,
        """
        Gravity config:
               Gravity Model: $(traitstring(config.model))
            Force Softenings: 
                         GAS => $(config.ForceSofteningTable[Int(GAS)])
                        STAR => $(config.ForceSofteningTable[Int(STAR)])
                        HALO => $(config.ForceSofteningTable[Int(HALO)])
                        DISL => $(config.ForceSofteningTable[Int(DISK)])
                       BULGE => $(config.ForceSofteningTable[Int(BULGE)])
                   BLACKHOLE => $(config.ForceSofteningTable[Int(BLACKHOLE)])


             OutboundLimiter: $(traitstring(config.outbound))
                 EnlargeMesh: $(config.enlarge)
           BoundaryCondition: $(traitstring(config.boundary))
                Sparse array: $(config.sparse)
            Index of Mondian interpolation function (`MOND_nuIndex`): $(config.MOND_nuIndex)
        """
    )
end


##### SolverConfig #####

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct SolverConfig{SolverG #=, SolverH=#}
    "Gravitational force solver. Supported: `DirectSum`, `Tree`, `FDM`, `FFT`, `ML`. Default is `DirectSum`"
    grav::SolverG
    #"Hydrodynamical force solver. Supported: `SPH`, `MHD`, `FEM`, `FVM`"
    #hydro::SolverH
end

"""
$(TYPEDSIGNATURES)
"""
function SolverConfig(;
    # Keywords from SimConfig
    GravitySolver::Gravity = DirectSum(),

    # Keywords to override
    grav::Gravity = GravitySolver,
)
    return SolverConfig(
        grav,
    )
end

function Base.show(io::IO, config::SolverConfig)
    print(io,
        """
        Solver config:
            Gravity solver: $(traitstring(config.grav))
        """
    )
end



##### DeviceConfig #####

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct DeviceConfig{Device, GA}
    "Supported: `CPU`, `GPU`"
    type::Device

    # for GPU
    GPU_NumThreads::Int
    "`::GPUAlgorithm`. Supported: `AllPairs`, `Tiled`. Default: `AllPairs`"
    GPU_Algorithm::GA
end

"""
$(TYPEDSIGNATURES)
"""
function DeviceConfig(;
    type::DeviceType = CPU(),
    GPU_NumThreads::Int = 32,
    GPU_Algorithm::GPUAlgorithm = AllPairs(),
)
    return DeviceConfig(
        type,
        GPU_NumThreads,
        GPU_Algorithm,
    )
end

function Base.show(io::IO, config::DeviceConfig)
    if config.type isa CPU
        print(io,
            """
            Device config:
                Type: CPU
                Number of threads: $(Threads.nthreads())
            """
        )
    elseif config.type isa GPU
        print(io,
            """
            Device config:
                             Type: GPU
                Number of threads: $(config.GPU_NumThreads)
            """
        )
    end
end



##### TimeConfig #####
abstract type AbstractTimestepConfig{T} end

struct AdaptiveTimestep{T<:Number,F} <: AbstractTimestepConfig{T}
    MinStep::T
    MaxStep::T
    TimeBase::Int64
    TimeInterval::T
    ErrTolTimestep::F
end

function Base.show(io::IO, config::AdaptiveTimestep)
    print(io,
        """
        AdaptiveTimestep:
            MinStep        = $(config.MinStep)
            MaxStep        = $(config.MaxStep)
            TimeBase       = $(config.TimeBase)
            TimeInterval   = $(config.TimeInterval)
            ErrTolTimestep = $(config.ErrTolTimestep)
        """
    )
end

struct ConstantTimestep{T<:Number} <: AbstractTimestepConfig{T}
    dt::T
end

function Base.show(io::IO, config::ConstantTimestep)
    print(io,
        """
        Constant Timestep:
            dt = $(config.dt)
        """
    )
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

## Notes on Keywords

- Keywords are compatible with `SimConfig`, and can also be overrided by field names of TimeConfig
- If timestep is non-zero, use `ConstantTimestep`; otherwise `AdaptiveTimesteps`
- `TimeBase` is the length of integer timeline. Default is `1<<30`
- `ErrTolTimestep` controls the accuracy of adaptive time integration

## Examples

```julia
TimeConfig(; TimeEnd = 1.0u"Gyr") # Keyword `End = TimeEnd` by default
TimeConfig(; End = 2.0u"Gyr")     # Override keyword from `SimConfig`
```
"""
struct TimeConfig{Time, Step, I, TIA}
    "Physical time at the beginning of the simulation"
    Begin::Time
    "Physical time at the end of the simulation"
    End::Time
    
    "Time interval between snapshots"
    BetweenSnapshots::Time
    "Integer time interval between snapshots"
    BetweenSnapshotsInt::I
    "Timestep config. Controled by keyword `TimeStep`. If `TimeStep` is zero, use `AdaptiveTimestep`; otherwise use `ConstantTimestep`"
    step::Step
    "Time integration algorithm. Supported: `Euler`, `Leapfrog`"
    algorithm::TIA
    "Redshift at start"
    redshift::Float64
    "Scale factor at start"
    scalefactor::Float64
end

"""
$(TYPEDSIGNATURES)
"""
function TimeConfig(;
    # Keywords from SimConfig
    TimeBegin::Number = 0.0u"Gyr",
    TimeEnd::Number = 0.1u"Gyr",
    TimeBetweenSnapshots::Number = 0.01u"Gyr",
    TimeStep::Number = 0.0u"Gyr",
    TimeBase = 1<<30,
    ErrTolTimestep = 0.025,
    MinStep::Number = isnothing(units) ? 1.0e-8 : 1.0e-8 * units[2],
    MaxStep::Number = isnothing(units) ? Inf : Inf * units[2],
    TimeIntegrationAlgorithm = Leapfrog(),

    redshift::Float64 = 0.0,
    scalefactor::Float64 = 1.0,

    # Keywords to override
    Begin::Number = TimeBegin,
    End::Number = TimeEnd,
    BetweenSnapshots::Number = TimeBetweenSnapshots,
    Step::Number = TimeStep,
    algorithm::TimeIntegration = TimeIntegrationAlgorithm,
)
    if iszero(ustrip(Step))
        if TimeBase < 1<<10
            @warn "TimeBase is too small!"
        end
        TimeInterval = (End - Begin) / TimeBase
        TimestepConfig = AdaptiveTimestep(MinStep, MaxStep, TimeBase, TimeInterval, ErrTolTimestep)
        TimeBetweenSnapshotsInt = floor(Int64, BetweenSnapshots / TimeInterval)
    else
        TimestepConfig = ConstantTimestep(Step)
        TimeBetweenSnapshotsInt = 0
    end

    return TimeConfig(
        Begin, End,
        BetweenSnapshots, TimeBetweenSnapshotsInt,
        TimestepConfig,
        algorithm,
        redshift, scalefactor,
    )
end

function Base.show(io::IO, config::TimeConfig)
    print(io,
        """
        Time config:
                             Begin: $(config.Begin)
                               End: $(config.End)
            Time between snapshots: $(config.BetweenSnapshots)
             Integration algorithm: $(traitstring(config.algorithm))
        """,
        tabstring(string(config.step)),
    )
end



##### SimConfig #####

"""
$(TYPEDEF)

$(TYPEDFIELDS)

## Examples
```julia
SimConfig(; GravitySolver = Tree())
SimConfig(; device = GPU(), TimeStep = 1.0e-5u"Gyr")
SimConfig(; units = uGadget2)
SimConfig(; units = nothing, ForceSofteningTable = [0.01 for i in 1:6])
SimConfig(; TimeEnd = 1.0u"Gyr", OutputDir = "Test/Dir")
```
"""
struct SimConfig{F, U}
    name::String
    author::String
    daytime::DateTime

    "Numeric type of float numbers"
    floattype::F

    "Prefered units in simulation. Default is `uAstro`. See `PhysicalParticles`: `uAstro`, `uSI`, `uGadget2`, `uCGS`. To run without units, use `nothing`"
    units::U
    "`ZeroValue`. Pre-constructed zero values for different types to simplify function arguments"
    ZeroValues
    "`Constant`. Physical constants"
    constants

    "Choose how to display status of simulation. Supported: `NormalMode`, `ProgressMode`, `SilentMode`. Default is `ProgressMode` displaying progress bars"
    loggingmode

    "`TimeConfig`"
    time
    "`OutputConfig`"
    output
    "`SolverConfig`"
    solver
    "`GravityConfig`"
    grav
    "`DeviceConfig`"
    device
end

"""
$(TYPEDSIGNATURES)
"""
function SimConfig( ;
    name = "Test",
    author = "AstroNbodySim",
    daytime = now(),
    floattype = Float64,
    units = uAstro,
    ZeroValues = ZeroValue(units),
    constants = Constant(units),
    loggingmode = ProgressMode(),

    # time
    TimeBegin = isnothing(units) ? floattype(0.0) : floattype(0.0) * units[2],
    TimeEnd = isnothing(units) ? floattype(0.1) : floattype(0.1) * units[2],
    TimeBetweenSnapshots = isnothing(units) ? floattype(0.001) : floattype(0.001) * units[2],
    TimeStep = isnothing(units) ? floattype(0.0) : floattype(0.0) * units[2],
    TimeBase = 1<<30,
    ErrTolTimestep = floattype(0.025),
    MinStep = isnothing(units) ? floattype(1.0e-8) : floattype(1.0e-8) * units[2],
    MaxStep = isnothing(units) ? typemax(floattype) : typemax(floattype) * units[2],
    TimeIntegrationAlgorithm = Leapfrog(),

    # output
    OutputDir = joinpath(pwd(), "output"),
    OutputFunction = output,
    OutputType = gadget2(),
    format2 = true,
    acc = false,
    pot = false,
    SaveRestart = false,
    SaveRestartFreq = 1.0,

    # solver
    GravitySolver = DirectSum(),

    # grav
    ForceSofteningTable = MVector{6}([
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
        isnothing(units) ? 1.0e-4 : 1.0e-4 * units[1],
    ]),
    MOND_nuIndex::Float64 = 2.0,
    GravityModel = Newton(),

    # device
    device::DeviceType = CPU(),
    GPU_NumThreads::Int = 32,
    GPU_Algorithm::GPUAlgorithm = AllPairs(),

    # mesh
    OutboundLimiter::OutboundLimiter = DS(),
    EnlargeMesh::Float64 = 2.01,
    BoundaryCondition::BoundaryCondition = Vacuum(),
    sparse::Bool = true,
)
    
    return SimConfig(
        name, author, daytime,
        floattype,
        units, ZeroValues, constants,
        loggingmode,
        TimeConfig(;
            TimeBegin, TimeEnd, TimeBetweenSnapshots, TimeStep, TimeBase,
            ErrTolTimestep, MinStep, MaxStep, TimeIntegrationAlgorithm,
        ),
        OutputConfig(;
            OutputDir, OutputFunction, OutputType,
            format2, acc, pot,
            SaveRestart, SaveRestartFreq,
        ),
        SolverConfig(;
            GravitySolver,
        ),
        GravityConfig(;
            ForceSofteningTable = MVector{length(ForceSofteningTable)}(ForceSofteningTable),
            GravityModel,
            MOND_nuIndex,
            OutboundLimiter,
            EnlargeMesh,
            BoundaryCondition,
            sparse,
        ),
        DeviceConfig(;
            type = device,
            GPU_NumThreads,
            GPU_Algorithm,
        ),
    )
end

function Base.show(io::IO, config::SimConfig)
    print(io,
        """
        Simulation Config: 
        """,
        tabstring(string(config.device)), "\n",
        tabstring(string(config.solver)), "\n",
        tabstring(string(config.grav)), "\n",
        tabstring(string(config.time)), "\n",
        tabstring(string(config.output)), "\n",
    )
end



##### LogInfo #####
@enum DefaultTimer begin
       TOTAL = 1
       FORCE = 2
       DRIFT = 3
        KICK = 4
    ANALYSIS = 5
      OUTPUT = 6
        PLOT = 7
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

## Examples #TODO
Use uppercase letters to avoid
```julia

```
"""
struct LogInfo
    "Timer enum names to access timing (continuously starting from 1)"
    timers::Symbol
    "Timings defined by timer enums"
    timing::Vector{UInt64}
    "Analyse on the whole simulation"
    analysers::Dict{String,Function}
end

"""
$(TYPEDSIGNATURES)
"""
function LogInfo(;
        timers = :DefaultTimer,
        analysers = Dict{String, Function}(),
    )
    timing = [UInt64(0) for i in instances(eval(timers))]
    return LogInfo(timers, timing, analysers)
end

"""
    Base.setproperty!(x::LogInfo, symbol::Symbol, d::Dict)

Modify the memories of Dict rather than modifing the pointer of Dict
"""
function Base.setproperty!(x::LogInfo, symbol::Symbol, d::Dict)
    a = getproperty(x, symbol)
    if !(a===d)
        empty!(a)
        merge!(a, d)
    end
end


##### StreamInfo #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct StreamInfo
    loggingio::IOStream
    timerio::IOStream
    analyserio::IOStream
end

StreamInfo() = StreamInfo(IOStream(""), IOStream(""), IOStream(""))


##### TimeInfo #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct TimeInfo{T<:Number, I<:Integer}
    "Time interval between neighbor time steps"
    dt::T

    system_time_int::I
    system_time_float::T
    last_system_time_int::I
    last_system_time_float::T
    
    next_output_time_int::I
    next_output_time_float::T

    redshift::Float64
    scalefactor::Float64

    min_endstep::I

    stepcount::I
end

"""
$(TYPEDSIGNATURES)
"""
function TimeInfo(config::SimConfig)
    dt = (isnothing(config.units) ? config.floattype(0.0) : config.floattype(0.0) * config.units[2])

    next_output_time_int = config.time.BetweenSnapshotsInt
    next_output_time_float = config.time.Begin + config.time.BetweenSnapshots

    system_time_float = config.time.Begin
    last_system_time_int = 0
    last_system_time_float = isnothing(config.units) ? config.floattype(0.0) : config.floattype(0.0) * config.units[2]

    return TimeInfo(
        dt,
        0, system_time_float,
        last_system_time_int, last_system_time_float,
        next_output_time_int, next_output_time_float,
        config.time.redshift,
        config.time.scalefactor,
        0,
        0,
    )
end

function Base.show(io::IO, info::TimeInfo)
    print(io,
    """
    Time Info:
                              dt: $(info.dt)
             system time Integer: $(info.system_time_int)
             system time Float  : $(info.system_time_float)
        last system time Integer: $(info.last_system_time_int)
        last system time Float  : $(info.last_system_time_float)
        next output time Integer: $(info.next_output_time_int)
        next output Time Float  : $(info.next_output_time_float)
                        redshift: $(info.redshift)
                     scalefactor: $(info.scalefactor)
                minimum end_step: $(info.min_endstep)
                       stepcount: $(info.stepcount)
    """
    )
end


##### OutputInfo #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct OutputInfo
    snapshotcount::Int64

    syncflag::Bool
end

"""
$(TYPEDSIGNATURES)
"""
function OutputInfo(;
        snapshotcount = 0,
    )
    return OutputInfo(
        snapshotcount,
        false,
    )
end

function Base.show(io::IO, info::OutputInfo)
    print(io,
    """
    Output Info:
        snapshotcount: $(info.snapshotcount)
             syncflag: $(info.syncflag)
    """
    )
end


##### PhysicsInfo #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct PhysicsInfo
    NumForceUpdateSinceLast::Int64
end

"""
$(TYPEDSIGNATURES)
"""
function PhysicsInfo()
    return PhysicsInfo(
        0,
    )
end

function Base.show(io::IO, info::PhysicsInfo)
    print(io,
    """
    Physics Info:
        NumForceUpdateSinceLast: $(info.NumForceUpdateSinceLast)
    """
    )
end



##### VisualizationInfo #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct VisualizationInfo
    progress::Progress

    PlotData
    resolution
    fig

    Realtime::Bool
    RenderTime::Float64   # s. Interactively render a new frame for every RenderTime
    last_plot_time::Float64 # s.

    xlims
    ylims
    zlims

    markersize::Float64
end

"""
$(TYPEDSIGNATURES)
"""
function VisualizationInfo(;
        Realtime::Bool = false,
        RenderTime::Float64 = 0.2,
        resolution = (1000, 1000),
        xlims = nothing,
        ylims = nothing,
        zlims = nothing,
        markersize = 0.0,
    )
    return VisualizationInfo(
        Progress(0),

        nothing,
        resolution,
        nothing,

        Realtime, RenderTime, 0.0,

        xlims, ylims, zlims, markersize,
    )
end



##### Simulation #####
"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Simulation{D}
    config
    
    id::Pair{Int,Int}
    pids::Array{Int,1}
    simdata::D

    # info
    "`::TimeInfo`"
    timeinfo
    "`::OutputInfo`"
    outputinfo
    "`::LogInfo`"
    loginfo
    "`::PhysicsInfo`"
    physics
    "`::StreamInfo`"
    stream
    "`::VisualizationInfo`"
    visinfo

    "`::Buffer`"
    buffer

    bgforce::Vector{Function}
    bgpotential::Vector{Function}
end

function Base.show(io::IO, sim::Simulation)
    if sim.config.solver.grav isa DirectSum && sim.config.device.type isa GPU
        DataCuts = numlocal(sim)
    else
        DataCuts = gather(sim, numlocal)
    end

    print(io,
        """
        Simulation Config:
                   id: $(sim.id)
                 pids: $(sim.pids)
            data cuts: $(DataCuts)
        """, "\n",
        tabstring(string(sim.config)), "\n",
        tabstring(string(sim.timeinfo)), "\n",
        tabstring(string(sim.outputinfo)), "\n",
        tabstring(string(sim.physics)), "\n",
        """
            Background force fields: $(sim.bgforce)

            Background potential fields: $(sim.bgpotential)
        """,
    )
end

"""
$(TYPEDSIGNATURES)
"""
function Simulation(d;
    floattype = Float64,
    units = uAstro,
    pids = [1],

    # LogInfo
    timers = :DefaultTimer,
    analysers = Dict{String, Function}(),

    # VisualizationInfo
    Realtime::Bool = false,
    RenderTime::Float64 = 0.2,
    resolution = (1000, 1000),

    xlims = nothing,
    ylims = nothing,
    zlims = nothing,
    markersize = 0.0,

    # background fields
    bgforce = Function[],
    bgpotential = Function[],

    # Tree method
    TreeOpenAngle = 0.1,
    ErrTolAcc = floattype(0.025),

    ToptreeAllocFactor::Int64 = 1,
    MaxTopnode = 20000,
    TopnodeFactor = 20,

    TreeAllocFactor::Int64 = 1,
    MaxTreenode = length(d) > 1000 ? length(d) * 10 : 10000,

    ExtentMargin = 1.001,

    PeanoBits3D = 21,
    PeanoBits2D = 31,

    epsilon = 1.0e-4,

    # Mesh setup
    meshmode = VertexMode(),
    assignment = CIC(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NG = 1,
    xMin = nothing,
    xMax = nothing,
    yMin = nothing,
    yMax = nothing,
    zMin = nothing,
    zMax = nothing,
    device = CPU(),
    EnlargeMesh = 2.01,
    BoundaryCondition = Vacuum(),

    # Other keywords are passed to SimConfig
    kw...
)
    config = SimConfig(; floattype, units, device, EnlargeMesh, BoundaryCondition, kw...)
    check_compatibility(config)
    id = next_simulationid()
    
    if config.solver.grav isa DirectSum
        if config.device.type isa CPU
            @info "Setting up Direct Sum simulation..."
            dStruct = StructArray(d)
            @sync for i in 1:length(pids)
                simdata = split_data(dStruct, i, length(pids))
                Distributed.remotecall_eval(AstroNbodySim, pids[i], :(
                    registry[$id] = Simulation(
                        $config, $id, $pids, $simdata,
                        TimeInfo($config),
                        OutputInfo(),
                        LogInfo(; timers = Symbol($timers), analysers = $analysers),
                        PhysicsInfo(),
                        StreamInfo(),
                        VisualizationInfo(;
                            Realtime = $Realtime, RenderTime = $RenderTime, resolution = $resolution,
                            xlims = $(xlims), ylims = $(ylims), zlims = $(zlims), markersize = $(markersize),
                        ),
                        Buffer($config),
                        $bgforce, $bgpotential,
                    )))
            end
            if !haskey(registry, id) # Master is not in pids, init empty data
                registry[id] = Simulation(
                    config, id, pids, 
                    empty(dStruct),
                    TimeInfo(config),
                    OutputInfo(),
                    LogInfo(; timers, analysers),
                    PhysicsInfo(),
                    StreamInfo(),
                    VisualizationInfo(; Realtime, RenderTime, resolution, xlims, ylims, zlims, markersize),
                    Buffer(config),
                    bgforce, bgpotential,
                )
            end
            @info "Data cuts: " * string(gather(registry[id], numlocal))
        elseif config.device.type isa GPU
            @info "Setting up Direct Sum simulation on GPU..."
            if config.solver.grav isa DirectSum
                dStruct = StructArray(d)
                registry[id] = Simulation(
                    config, id, pids,
                    cu(dStruct),
                    TimeInfo(config),
                    OutputInfo(),
                    LogInfo(; timers, analysers),
                    PhysicsInfo(),
                    StreamInfo(),
                    VisualizationInfo(; Realtime, RenderTime, resolution, xlims, ylims, zlims, markersize),
                    Buffer(config),
                    bgforce, bgpotential,
                )
            end
        end
    elseif config.solver.grav isa Tree
        @info "Setting up Tree simulation..."
        octreeconfig = OctreeConfig(length(d);
            ToptreeAllocFactor,
            MaxTopnode,
            TopnodeFactor,
            TreeAllocFactor,
            MaxTreenode,
            ExtentMargin,
            PeanoBits3D,
            PeanoBits2D,
            epsilon,
        )
        tree = octree(d, config = octreeconfig, pids = pids, units = units)
        
        treesimconfig = TreeSimConfig(;
            TreeOpenAngle,
            ErrTolAcc,
        )

        @everywhere pids AstroNbodySim.registry[$id] = Simulation(
            $config, $id, $pids,
            OctreeData($treesimconfig, AstroNbodySim.PhysicalTrees.registry[$(tree.id)]),
            TimeInfo($config),
            OutputInfo(),
            LogInfo(; timers = Symbol($timers), analysers = $analysers),
            PhysicsInfo(),
            StreamInfo(),
            VisualizationInfo(;
                Realtime = $Realtime, RenderTime = $RenderTime, resolution = $resolution,
                xlims = $(xlims), ylims = $(ylims), zlims = $(zlims), markersize = $(markersize),
            ),
            Buffer($config),
            $bgforce, $bgpotential,
        )

        if !haskey(registry, id) # Master is not in pids, init empty data
            registry[id] = Simulation(
                config, id, pids,
                OctreeData(treesimconfig, tree),
                TimeInfo(config),
                OutputInfo(),
                LogInfo(; timers, analysers),
                PhysicsInfo(),
                StreamInfo(),
                VisualizationInfo(; Realtime, RenderTime, resolution, xlims, ylims, zlims, markersize),
                Buffer(config),
                bgforce, bgpotential,
            )
        end
        @info "Data cuts: " * string(gather(registry[id], numlocal))
    elseif config.solver.grav isa FDM || config.solver.grav isa FFT
        @info "Setting up $(traitstring(config.solver.grav)) simulation..."
        dStruct = StructArray(d)
        mesh = MeshCartesianStatic(dStruct, units;
            boundary = BoundaryCondition,
            assignment,
            Nx, Ny, Nz, NG,
            xMin, xMax, yMin, yMax, zMin, zMax,
            mode = meshmode,
            gpu = device isa GPU ? true : false,
            enlarge = EnlargeMesh,
        )
        registry[id] = Simulation(
            config, id, pids,
            mesh,
            TimeInfo(config),
            OutputInfo(),
            LogInfo(; timers, analysers),
            PhysicsInfo(),
            StreamInfo(),
            VisualizationInfo(; Realtime, RenderTime, resolution, xlims, ylims, zlims, markersize),
            Buffer(config),
            bgforce, bgpotential,
        )
    elseif config.solver.grav isa ML
    end

    return registry[id]
end

@inline Base.length(::Simulation) = 1
@inline Base.iterate(s::Simulation) = (s, nothing)
@inline Base.iterate(s::Simulation, st) = nothing

function getunits(sim::Simulation)
    return sim.config.units
end

function get_local_data(sim::Simulation)
    return get_local_data(sim, sim.config.solver.grav, sim.config.device.type)
end

function get_local_data(sim::Simulation, ::DirectSum, ::CPU)
    return sim.simdata
end

function get_local_data(sim::Simulation, ::Tree, ::CPU)
    return sim.simdata.tree.data
end

function get_local_data(sim::Simulation, ::Union{FDM, FFT}, ::DeviceType)
    return sim.simdata.data
end

"Access data on GPU"
function get_local_data(sim::Simulation, ::DirectSum, ::GPU)
    return sim.simdata
end


function numlocal(sim::Simulation)
    return length(get_local_data(sim))
end


function get_all_data(sim::Simulation)
    return get_all_data(sim, sim.config.solver.grav, sim.config.device.type)
end

function get_all_data(sim::Simulation, ::Union{DirectSum, Tree}, ::CPU)
    data = gather(sim, get_local_data)
    d = deepcopy(first(data))

    for p in data[2:end]
        append!(d, p)
    end
    return d
end

function get_all_data(sim::Simulation, ::Union{FDM, FFT}, ::CPU)
    return sim.simdata.data
end

"Copy data to CPU"
function get_all_data(sim::Simulation, ::DirectSum, ::GPU)
    CUDA.@allowscalar return StructArray(Array(sim.simdata))
end

"Copy data to CPU"
function get_all_data(sim::Simulation, ::Union{FDM, FFT}, ::GPU)
    CUDA.@allowscalar return StructArray(Array(sim.simdata.data))
end

function check_QUMOND(config::SimConfig)
    incompatible = config.grav.model isa QUMOND && config.solver.grav isa FDM && !(config.grav.boundary isa Vacuum)
    if incompatible
        error("QUMOND can only be solved by FDM method in the vacuum boundary condition")
    end
end

function check_tree(config::SimConfig)
    incompatible = config.solver.grav isa Tree && config.device.type isa GPU
    if incompatible
        error("Tree simulation on GPU is not supported!")
    end
end

function check_mesh(config::SimConfig)
    incompatible = config.solver.grav isa FDM && config.device.type isa GPU && config.grav.sparse
    if incompatible
        error("Cannot use SparseArray on GPU in FDM simulation")
    end

    incompatible = config.solver.grav isa Union{FDM,FFT} && config.device.type isa GPU
    if incompatible
        @warn "Dynamic Particle-mesh simulation is not suggested to use GPU. It is still in DEV stage."
    end
end

function check_compatibility(config::SimConfig)
    check_QUMOND(config)
    check_tree(config)
    check_mesh(config)
end