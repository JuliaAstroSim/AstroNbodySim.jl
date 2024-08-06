mkpathIfNotExist("output/energy")

# Load or generate data
if isfile("output/energy/EnergyMomentum.gadget2")
    header, _d = read_gadget2("output/energy/EnergyMomentum.gadget2", uAstro, type=Star);
    _d64 = convert.(Float64, _d)
    particles = convert.(Int, _d64)
else
    config = PlummerStarCluster(
        collection = STAR,
        NumSamples = 500,
        VirialRadius = 0.010u"kpc",
        TotalMass = 1.0e5u"Msun",
        model = AstroIC.Newton(),
    )
    particles = generate(config, MaxRadius = 0.2u"kpc");
    write_gadget2("output/energy/EnergyMomentum.gadget2", particles, format2 = false)
end

TimeEnd = 0.05u"Gyr"
TimeBetweenSnapshots = 0.005u"Gyr"

analysers = Dict{String, Function}(
    "potential" => total_potential,
    "kinetic" => total_kinetic,
    "momentum" => total_momentum,
    "angularmomentum" => total_angular_momentum,
)


# The conservation of momentum is dependent upon time integration accuracy
# There are two ways to achieve a better level of conservation
#   1. increase SoftLen to suppress close collisions and binaries
#   2. decrease ErrTolTimestep for higher time integration accuracy (ErrTolIntAccuracy in Gadget2 param)

#=
SoftLen = 0.0001u"kpc"
ForceSofteningTable = [SoftLen for i in 1:6]
ErrTolTimestep = 0.025
=#

function test_energy_momentum(SoftLen::Number, ErrTolTimestep::Float64, TimeStep::Number;
    gadget = false,
    plot_only = false,
)
    ForceSofteningTable = [SoftLen for i in 1:6]
    
    @info "Setting up simulations"
    ds = Simulation(
        deepcopy(particles);
        TimeEnd,
        TimeBetweenSnapshots,
        TimeStep,
        analysers,
        ForceSofteningTable,
        ErrTolTimestep,
        OutputDir = "output/energy/DirectSum_$(ustrip(SoftLen))_$(ErrTolTimestep)",
    );

    ds_adapt = Simulation(
        deepcopy(particles);
        TimeEnd,
        TimeBetweenSnapshots,
        TimeStep = 0.0u"Gyr",
        analysers,
        ForceSofteningTable,
        ErrTolTimestep,
        OutputDir = "output/energy/DirectSumAdaptive_$(ustrip(SoftLen))_$(ErrTolTimestep)",
    );

    ts = Simulation(
        deepcopy(particles);
        TimeEnd,
        TimeBetweenSnapshots,
        TimeStep,
        analysers,
        ForceSofteningTable,
        ErrTolTimestep,
        GravitySolver = Tree(),
        OutputDir = "output/energy/Tree_$(ustrip(SoftLen))_$(ErrTolTimestep)",
    );

    ts_adapt = Simulation(
        deepcopy(particles);
        TimeEnd,    
        TimeBetweenSnapshots,
        TimeStep = 0.0u"Gyr",
        analysers,
        ForceSofteningTable,
        ErrTolTimestep,
        GravitySolver = Tree(),
        OutputDir = "output/energy/TreeAdaptive_$(ustrip(SoftLen))_$(ErrTolTimestep)",
    );

    if !plot_only
        @info "Running simulations"
        run(ds)
        run(ds_adapt)
        run(ts)
        run(ts_adapt)
    end

    @info "Plotting results"
    fig1, df1 = plot_energy(joinpath(ds.config.output.dir, "analysis.csv"), size=(800,450), margin = (10, 10, 10, 50));
    fig2, df2 = plot_energy_delta(joinpath(ds.config.output.dir, "analysis.csv"), size=(800,450));
    Makie.save("output/energy/DirectSum_$(ustrip(SoftLen))_$(ErrTolTimestep)_Energy.png", fig1)
    Makie.save("output/energy/DirectSum_$(ustrip(SoftLen))_$(ErrTolTimestep)_DeltaEnergy.png", fig2)

    fig3, df3 = plot_energy(joinpath(ds_adapt.config.output.dir, "analysis.csv"), size=(800,450), margin = (10, 10, 10, 50));
    fig4, df4 = plot_energy_delta(joinpath(ds_adapt.config.output.dir, "analysis.csv"), size=(800,450));
    Makie.save("output/energy/DirectSumAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_Energy.png", fig3)
    Makie.save("output/energy/DirectSumAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_DeltaEnergy.png", fig4)

    fig5, df5 = plot_energy(joinpath(ts.config.output.dir, "analysis.csv"), size=(800,450), margin = (10, 10, 10, 50));
    fig6, df6 = plot_energy_delta(joinpath(ts.config.output.dir, "analysis.csv"), size=(800,450));
    Makie.save("output/energy/Tree_$(ustrip(SoftLen))_$(ErrTolTimestep)_Energy.png", fig5)
    Makie.save("output/energy/Tree_$(ustrip(SoftLen))_$(ErrTolTimestep)_DeltaEnergy.png", fig6)

    fig7, df7 = plot_energy(joinpath(ts_adapt.config.output.dir, "analysis.csv"), size=(800,450), margin = (10, 10, 10, 50));
    fig8, df8 = plot_energy_delta(joinpath(ts_adapt.config.output.dir, "analysis.csv"), size=(800,450));
    Makie.save("output/energy/TreeAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_Energy.png", fig7)
    Makie.save("output/energy/TreeAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_DeltaEnergy.png", fig8)


    # For constant timesteps, they are having `range step cannot be zero` error from Makie, try Plots.jl instead
    # The reason is that they are conserved at smallest 32-bit floating point level. So we omit the plotting.
    # Ref: https://github.com/JuliaPlots/Makie.jl/issues/1579
    #fig1, df1 = plot_momentum(joinpath(ds.config.output.dir, "analysis.csv"); colors, size=(800,450));
    #fig2, df2 = plot_momentum_angular(joinpath(ds.config.output.dir, "analysis.csv"); colors, absolute = true, size=(800,450));
    #Makie.save("output/energy/DirectSum_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.png", fig1)
    #Makie.save("output/energy/DirectSum_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum_angular.png", fig2)

    fig3, df3 = plot_momentum(joinpath(ds_adapt.config.output.dir, "analysis.csv"); colors, size=(800,450));
    fig4, df4 = plot_momentum_angular(joinpath(ds_adapt.config.output.dir, "analysis.csv"); colors, absolute = true, size=(800,450));
    Makie.save("output/energy/DirectSumAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.png", fig3)
    Makie.save("output/energy/DirectSumAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum_angular.png", fig4)

    #fig5, df5 = plot_momentum(joinpath(ts.config.output.dir, "analysis.csv"); colors, size=(800,450));
    #fig6, df6 = plot_momentum_angular(joinpath(ts.config.output.dir, "analysis.csv"); colors, absolute = true, size=(800,450));
    #Makie.save("output/energy/Tree_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.png", fig5)
    #Makie.save("output/energy/Tree_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum_angular.png", fig6)

    fig7, df7 = plot_momentum(joinpath(ts_adapt.config.output.dir, "analysis.csv"); colors, size=(800,450));
    fig8, df8 = plot_momentum_angular(joinpath(ts_adapt.config.output.dir, "analysis.csv"); colors, absolute = true, size=(800,450));
    Makie.save("output/energy/TreeAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.png", fig7)
    Makie.save("output/energy/TreeAdapt_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum_angular.png", fig8)


    # Compare Gadget
    if gadget
        @assert Base.Sys.isunix() "Gadget-2 test can only run in unix-like environments"
        if !isfile("output/energy/Gadget2/Gadget2/Gadget2")
            @info "Setting up Gadget2"
            rm("output/energy/Gadget2", force = true, recursive = true)
            mkdir("output/energy/Gadget2")
            cd("output/energy/Gadget2")
            
            # compile Gadget2
            cp(joinpath(@__DIR__, "../Gadget2"), "Gadget2", force = true)
            cd("Gadget2")
            write_gadget2_makefile("Makefile", "OUTPUTPOTENTIAL")
            
            # Make sure that gravity constant is the same
            G_cgs = ustrip(uconvert(u"cm^3/g/s^2", Constant().G))
            alter_line("allvars.h", "#define  GRAVITY",
            "#define  GRAVITY           $(G_cgs)   /*!< Gravitational constant (in cgs units) */")
            
            alter_line("io.c", """sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);""",
            """sprintf(buf, "%s%s_%04d", All.OutputDir, All.SnapshotFileBase, num);""");
            
            result = run(`make -j 4`)
            cd("../")
        else
            cd("output/energy/Gadget2")
        end
        
        rm("output_$(ustrip(SoftLen))_$(ErrTolTimestep)", force = true, recursive = true)
        mkdir("output_$(ustrip(SoftLen))_$(ErrTolTimestep)")
    
        # Output in format2, which is easier to read POT
        write_gadget2_param("IC_$(ustrip(SoftLen))_$(ErrTolTimestep).param", "../EnergyMomentum.gadget2", ts;
            OutputDir = "output_$(ustrip(SoftLen))_$(ErrTolTimestep)",
            SnapFormat = 1,
            TimeBetSnapshot = 0.0001,
            ErrTolIntAccuracy = ErrTolTimestep,
        )
    
        result = run(`mpirun -np 4 Gadget2/Gadget2 IC_$(ustrip(SoftLen))_$(ErrTolTimestep).param`)
        cd(@__DIR__)
    
        # Analyse data
    
        fig9, df9 = plot_energy(
            "output/energy/Gadget2/output_$(ustrip(SoftLen))_$(ErrTolTimestep)", "snapshot_", collect(0:500), "", gadget2(),
            times = collect(0.0u"Gyr":0.0001u"Gyr":TimeEnd),
            formatstring = "%04d",
            size=(800,450),
        )
        mv("energy.csv", "output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_energy.csv", force = true)
        Makie.save("output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_energy.png", fig9)
    
        fig10, df10 = plot_energy_delta("output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_energy.csv", size=(800,450))
        Makie.save("output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_energy_delta.png", fig10)
    
    
    
        fig11, df11 = plot_momentum(
            "output/energy/Gadget2/output_$(ustrip(SoftLen))_$(ErrTolTimestep)", "snapshot_", collect(0:500), "", gadget2(),
            times = collect(0.0u"Gyr":0.0001u"Gyr":TimeEnd),
            formatstring = "%04d",
            size=(800,450),
        )
        mv("momentum.csv", "output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.csv", force = true)
        Makie.save("output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum.png", fig11)
    
        fig12, df12 = plot_momentum_angular(df11; colors, absolute = true, size=(800,450))
        Makie.save("output/energy/Gadget2_$(ustrip(SoftLen))_$(ErrTolTimestep)_momentum_angular.png", fig12)
    end
    ISLENT.clear()
    return true
end

@testset "Energy & Momentum" begin
    @test test_energy_momentum(0.001u"kpc", 0.025, 0.000045u"Gyr")
    @test test_energy_momentum(0.01u"kpc", 0.025, 0.0003u"Gyr")
end

#=
test_energy_momentum(0.001u"kpc", 0.025, 4.5e-5u"Gyr")

test_energy_momentum(0.01u"kpc", 0.025, 0.0003u"Gyr")

test_energy_momentum(0.0001u"kpc", 0.025)
test_energy_momentum(0.001u"kpc", 0.025)
test_energy_momentum(0.01u"kpc", 0.025)


test_energy_momentum(0.0001u"kpc", 0.01)
test_energy_momentum(0.001u"kpc", 0.01)
test_energy_momentum(0.01u"kpc", 0.01)


test_energy_momentum(0.0001u"kpc", 0.005)
test_energy_momentum(0.001u"kpc", 0.005)
test_energy_momentum(0.01u"kpc", 0.005)



TimeStep = 4.5e-5u"Gyr"
test_energy_momentum(0.001u"kpc", 0.01)

TimeStep = 0.0001u"Gyr"
test_energy_momentum(0.01u"kpc", 0.01)
=#
