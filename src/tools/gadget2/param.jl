"""
function write_gadget2_param(filename::String, ICfilename::String)

# Params
- filename   : name of param file
- ICfilename : name of initial conditions
"""
function write_gadget2_param(filename::String, ICfilename::String;
    OutputDir = "output",
    SnapshotFileBase = "snapshot",
    TimeLimitCPU = 7200000,        # 7200000 s ~ 2000 hours

    # Code options
    ICFormat                = 1,  
    SnapFormat              = 1,
    ComovingIntegrationOn   = 0,
    TypeOfTimestepCriterion = 0,
    OutputListOn            = 0,
    PeriodicBoundariesOn    = 0,

    # Caracteristics of run
    TimeBegin   = 0.0,
    TimeMax     = 0.1,
    Omega0      = 0,
    OmegaLambda = 0,
    OmegaBaryon = 0,
    HubbleParam = 1.0,
    BoxSize     = 0,

    # Output frequency
    TimeBetSnapshot           = 0.0001,
    TimeOfFirstSnapshot       = 0,
    CpuTimeBetRestartFile     = 360000.0,
    TimeBetStatistics         = 0.05,
    NumFilesPerSnapshot       = 1,
    NumFilesWrittenInParallel = 1,

    # Accuracy of time integration
    ErrTolIntAccuracy      = 0.025,
    CourantFac             = 0.15,
    MaxSizeTimestep        = 0.01,
    MinSizeTimestep        = 0.000001,

    # Tree algorithm, force accuracy, domain update frequency
    ErrTolTheta                = 0.5,
    TypeOfOpeningCriterion     = 1,
    ErrTolForceAcc             = 0.005,

    # Further parameters of SPH
    DesNumNgb            = 64,
    MaxNumNgbDeviation   = 200,
    ArtBulkViscConst     = 0.8,
    InitGasTemp          = 0,
    MinGasTemp           = 0,

    # Memory allocation
    PartAllocFactor      = 3.5,
    TreeAllocFactor      = 2.8,
    BufferSize           = 100,

    # System of units (cgs)
    UnitLength_in_cm         = ustrip(uconvert(u"cm", 1.0u"kpc")),          #    1.0 kpc
    UnitMass_in_g            = 10^10 * ustrip(uconvert(u"g", 1.0u"Msun")),  # 1.0e10 Msun
    UnitVelocity_in_cm_per_s = 1.0e5,                                       #    1.0 km/s

    ForceSofteningTable    = [0.0001 for i in 1:6],
    ForceSofteningTableMax = [0.0001 for i in 1:6],
)
    gadgetparamstring = """
    InitCondFile  	   $ICfilename
    OutputDir          $OutputDir
    EnergyFile         energy.txt
    InfoFile           info.txt
    TimingsFile        timings.txt
    CpuFile            cpu.txt
    RestartFile        restart
    SnapshotFileBase   $SnapshotFileBase
    OutputListFilename parameterfiles/output_list.txt

    % CPU time -limit
    TimeLimitCPU      $TimeLimitCPU
    ResubmitOn        0
    ResubmitCommand   my-scriptfile  

    % Code options
    ICFormat                 $ICFormat
    SnapFormat               $SnapFormat
    ComovingIntegrationOn    $ComovingIntegrationOn
    TypeOfTimestepCriterion  $TypeOfTimestepCriterion
    OutputListOn             $OutputListOn
    PeriodicBoundariesOn     $PeriodicBoundariesOn
    
    %  Caracteristics of run
    TimeBegin           $TimeBegin
    TimeMax	            $TimeMax
    Omega0	            $Omega0
    OmegaLambda         $OmegaLambda
    OmegaBaryon         $OmegaBaryon
    HubbleParam         $HubbleParam
    BoxSize             $BoxSize
    
    % Output frequency
    TimeBetSnapshot           $TimeBetSnapshot
    TimeOfFirstSnapshot       $TimeOfFirstSnapshot
    CpuTimeBetRestartFile     $CpuTimeBetRestartFile
    TimeBetStatistics         $TimeBetStatistics
    NumFilesPerSnapshot       $NumFilesPerSnapshot
    NumFilesWrittenInParallel $NumFilesWrittenInParallel
    
    % Accuracy of time integration
    ErrTolIntAccuracy      $ErrTolIntAccuracy
    CourantFac             $CourantFac
    MaxSizeTimestep        $MaxSizeTimestep
    MinSizeTimestep        $MinSizeTimestep
    
    % Tree algorithm, force accuracy, domain update frequency
    ErrTolTheta                  $ErrTolTheta
    TypeOfOpeningCriterion       $TypeOfOpeningCriterion
    ErrTolForceAcc               $ErrTolForceAcc
    TreeDomainUpdateFrequency    0.0
    
    %  Further parameters of SPH
    DesNumNgb              $DesNumNgb
    MaxNumNgbDeviation     $MaxNumNgbDeviation
    ArtBulkViscConst       $ArtBulkViscConst
    InitGasTemp            $InitGasTemp
    MinGasTemp             $MinGasTemp
    
    % Memory allocation
    PartAllocFactor       $PartAllocFactor
    TreeAllocFactor       $TreeAllocFactor
    BufferSize            $BufferSize
    
    % System of units
    UnitLength_in_cm         $UnitLength_in_cm         #    1.0 kpc
    UnitMass_in_g            $UnitMass_in_g            # 1.0e10 Msun
    UnitVelocity_in_cm_per_s $UnitVelocity_in_cm_per_s #    1.0 km/s
    GravityConstantInternal  0

    % Softening lengths
    MinGasHsmlFractional 0.25
    SofteningGas       $(ForceSofteningTable[1])
    SofteningHalo      $(ForceSofteningTable[2])
    SofteningDisk      $(ForceSofteningTable[3])
    SofteningBulge     $(ForceSofteningTable[4])
    SofteningStars     $(ForceSofteningTable[5])
    SofteningBndry     $(ForceSofteningTable[6])

    SofteningGasMaxPhys       $(ForceSofteningTableMax[1])
    SofteningHaloMaxPhys      $(ForceSofteningTableMax[2])
    SofteningDiskMaxPhys      $(ForceSofteningTableMax[3])
    SofteningBulgeMaxPhys     $(ForceSofteningTableMax[4])
    SofteningStarsMaxPhys     $(ForceSofteningTableMax[5])
    SofteningBndryMaxPhys     $(ForceSofteningTableMax[6])
    
    MaxRMSDisplacementFac 0.2
    """
    f = open(filename, "w")
    write(f, gadgetparamstring)
    close(f)
    return true
end

"""
alter_param(filename::String, param::String, value)

alter all lines containing `param` to `param  value`
"""
function alter_param(filename::String, param::String, value)
    newfilename = filename * ".temp"
    newfile = open(newfilename, "w")

    words = readlines(filename)
    for line in words
        if occursin(param, line)
            write(newfile, "$param  $value")
        else
            write(newfile, line)
        end
        write(newfile, "\n")
    end

    close(newfile)
    mv(newfilename, filename, force = true)
end

"""
alter_line(filename::String, match::String, newline::String; all = false)

alter the first line containing `match` to `newline`

# Keywords
- all::Bool = false : If true, alter all lines containing `match` to `newline`
"""
function alter_line(filename::String, match::String, newline::String; all = false)
    newfilename = filename * ".temp"
    newfile = open(newfilename, "w")

    #TODO try blocks
    words = readlines(filename)
    firstdone = false
    for line in words
        if occursin(match, line)
            if !firstdone || all
                println(line)
                write(newfile, newline)
                firstdone = true
            end
        else
            write(newfile, line)
        end
        write(newfile, "\n")
    end

    close(newfile)
    mv(newfilename, filename, force = true)
end

function add_line(filename::String, newline::String)
    newfilename = filename * ".temp"
    newfile = open(newfilename, "w")

    words = readlines(filename)
    write(newfile, newline)
    write(newfile, "\n")
    for line in words
        write(newfile, line)
        write(newfile, "\n")
    end

    close(newfile)
    mv(newfilename, filename, force = true)
end

function write_gadget2_param(filename::String, ICfilename::String, config::SimConfig;
    kw...
)
    write_gadget2_param(filename, ICfilename,
        TimeBegin = ustrip(u"Gyr", config.time.Begin),
        TimeMax   = ustrip(u"Gyr", config.time.End),
        TimeBetSnapshot = ustrip(u"Gyr", config.time.BetweenSnapshots),
        ForceSofteningTable = ustrip.(u"kpc", config.grav.ForceSofteningTable),
        ; kw...
    )

    if config.time.step isa AdaptiveTimestep
        alter_param(filename, "ErrTolIntAccuracy", config.time.step.ErrTolTimestep)
    end
end

function write_gadget2_param(filename::String, ICfilename::String, sim::Simulation, ::Gravity, ::DeviceType; kw...)
    write_gadget2_param(filename, ICfilename, sim.config; kw...)
end

function write_gadget2_param(filename::String, ICfilename::String, sim::Simulation, ::Tree, ::CPU; kw...)
    write_gadget2_param(filename, ICfilename, sim.config;
        ErrTolTheta = sim.simdata.treesimconfig.TreeOpenAngle,
        ErrTolForceAcc = sim.simdata.treesimconfig.ErrTolAcc,
        kw...
    )
end

function write_gadget2_param(filename::String, ICfilename::String, sim::Simulation; kw...)
    write_gadget2_param(filename, ICfilename, sim, sim.config.solver.grav, sim.config.device.type; kw...)
end


function write_gadget4_param(filename::String, ICfilename::String;
    OutputDir = "output",
    SnapshotFileBase = "snapshot",
    TimeLimitCPU = 7200000,        # 7200000 s ~ 2000 hours

    # Code options
    ICFormat                = 1,  
    SnapFormat              = 1,
    ComovingIntegrationOn   = 0,
    TypeOfTimestepCriterion = 0,
    OutputListOn            = 0,

    # Caracteristics of run
    TimeBegin   = 0.0,
    TimeMax     = 0.1,
    Omega0      = 0,
    OmegaLambda = 0,
    OmegaBaryon = 0,
    HubbleParam = 1.0,
    BoxSize     = 0,

    # Output frequency
    TimeBetSnapshot           = 0.0001,
    TimeOfFirstSnapshot       = 0,
    CpuTimeBetRestartFile     = 360000.0,
    TimeBetStatistics         = 0.05,
    NumFilesPerSnapshot       = 1,

    # Accuracy of time integration
    ErrTolIntAccuracy      = 0.025,
    CourantFac             = 0.15,
    MaxSizeTimestep        = 0.01,
    MinSizeTimestep        = 0.000001,

    # Tree algorithm, force accuracy, domain update frequency
    ErrTolTheta                = 0.5,
    ErrTolForceAcc             = 0.005,

    # Further parameters of SPH
    DesNumNgb            = 64,
    MaxNumNgbDeviation   = 200,
    ArtBulkViscConst     = 0.8,
    InitGasTemp          = 0,
    MinGasTemp           = 0,

    # System of units (cgs)
    UnitLength_in_cm         = ustrip(uconvert(u"cm", 1.0u"kpc")),          #    1.0 kpc
    UnitMass_in_g            = 10^10 * ustrip(uconvert(u"g", 1.0u"Msun")),  # 1.0e10 Msun
    UnitVelocity_in_cm_per_s = 1.0e5,                                       #    1.0 km/s
)
    gadgetparamstring = """
    %----  Relevant files 
    InitCondFile  	   $ICfilename
    OutputDir          $OutputDir
    SnapshotFileBase    $SnapshotFileBase
    OutputListFilename  empty.txt

    %---- File formats
    ICFormat                 $ICFormat
    SnapFormat               $SnapFormat
    
    %---- CPU-time limits
    TimeLimitCPU             $TimeLimitCPU          % in seconds
    CpuTimeBetRestartFile    $CpuTimeBetRestartFile % in seconds

    %----- Memory alloction
    MaxMemSize        2300

    %---- Caracteristics of run
    TimeBegin           $TimeBegin
    TimeMax	            $TimeMax

    %---- Basic code options that set the type of sim
    ComovingIntegrationOn    $ComovingIntegrationOn
    
    %---- Cosmological parameters
    Omega0	            $Omega0
    OmegaLambda         $OmegaLambda
    OmegaBaryon         $OmegaBaryon
    HubbleParam         $HubbleParam
    Hubble              0.1
    BoxSize             $BoxSize
    
    %---- Output frequency and output paramaters
    OutputListOn             $OutputListOn
    TimeBetSnapshot           $TimeBetSnapshot
    TimeOfFirstSnapshot       $TimeOfFirstSnapshot
    TimeBetStatistics         $TimeBetStatistics
    NumFilesPerSnapshot       $NumFilesPerSnapshot
    MaxFilesWithConcurrentIO  1 

    % Accuracy of time integration
    ErrTolIntAccuracy      $ErrTolIntAccuracy
    CourantFac             $CourantFac
    MaxSizeTimestep        $MaxSizeTimestep
    MinSizeTimestep        $MinSizeTimestep

    %---- Tree algorithm, force accuracy, domain update frequency
    TypeOfTimestepCriterion  $TypeOfTimestepCriterion
    ErrTolTheta                  $ErrTolTheta
    ErrTolThetaMax                        1.0
    ErrTolForceAcc               $ErrTolForceAcc
    TopNodeFactor                         2.5

    ActivePartFracForNewDomainDecomp      0.01
    ActivePartFracForPMinsteadOfEwald     0.05

    %---- Initial density estimate
    DesNumNgb              $DesNumNgb
    MaxNumNgbDeviation     $MaxNumNgbDeviation

    %---- System of units
    UnitLength_in_cm         $UnitLength_in_cm         #    1.0 kpc
    UnitMass_in_g            $UnitMass_in_g            # 1.0e10 Msun
    UnitVelocity_in_cm_per_s $UnitVelocity_in_cm_per_s #    1.0 km/s
    GravityConstantInternal  0

    %---- Gravitational softening length
    SofteningComovingClass0      0
    SofteningComovingClass1      72.0
    SofteningComovingClass2      180.0
    SofteningComovingClass3      500.0
    
    SofteningMaxPhysClass0       0
    SofteningMaxPhysClass1       12.0
    SofteningMaxPhysClass2       30.0
    SofteningMaxPhysClass3       150.0
    
    SofteningClassOfPartType0    0
    SofteningClassOfPartType1    1
    SofteningClassOfPartType2    2
    SofteningClassOfPartType3    3
    
    %----- SPH
    ArtBulkViscConst       $ArtBulkViscConst
    InitGasTemp            $InitGasTemp
    MinGasTemp             $MinGasTemp
    """
    f = open(filename, "w")
    write(f, gadgetparamstring)
    close(f)
    return true
end

function write_gadget4_param(filename::String, ICfilename::String, config::SimConfig;
    kw...
)
    write_gadget4_param(filename, ICfilename,
        TimeBegin = ustrip(u"Gyr", config.time.Begin),
        TimeMax   = ustrip(u"Gyr", config.time.End),
        TimeBetSnapshot = ustrip(u"Gyr", config.time.BetweenSnapshots),
        ; kw...
    )

    if config.time.step isa AdaptiveTimestep
        alter_param(filename, "ErrTolIntAccuracy", config.time.step.ErrTolTimestep)
    end
end

function write_gadget4_param(filename::String, ICfilename::String, sim::Simulation, ::Gravity, ::DeviceType; kw...)
    write_gadget4_param(filename, ICfilename, sim.config; kw...)
end

function write_gadget4_param(filename::String, ICfilename::String, sim::Simulation, ::Tree, ::CPU; kw...)
    write_gadget4_param(filename, ICfilename, sim.config;
        ErrTolTheta = sim.simdata.treesimconfig.TreeOpenAngle,
        ErrTolForceAcc = sim.simdata.treesimconfig.ErrTolAcc,
        kw...
    )
end

function write_gadget4_param(filename::String, ICfilename::String, sim::Simulation; kw...)
    write_gadget4_param(filename, ICfilename, sim, sim.config.solver.grav, sim.config.device.type; kw...)
end