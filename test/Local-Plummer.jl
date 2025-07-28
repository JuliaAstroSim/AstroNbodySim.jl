#! This is a heavy-load test. Run at local computers!
mkpathIfNotExist("output/Plummer/")

function plot_plummer(sim::Simulation, mode)
    # Plot profiling
    fig = plot_profiling(joinpath(sim.config.output.dir, "timing.csv"),
        title = "Profiling of " * mode,
        size = (1600,900),
    )
    Makie.save("output/Plummer/Profiling-$(mode).png", fig)

    # Plot energy
    fig, df = plot_energy(joinpath(sim.config.output.dir, "analysis.csv"))
    Makie.save("output/Plummer/Energy-$(mode).png", fig)

    fig, df = plot_energy_delta(joinpath(sim.config.output.dir, "analysis.csv"))
    Makie.save("output/Plummer/EnergyDelta-$(mode).png", fig)

    # Plot momentum
    fig, df = plot_momentum(joinpath(sim.config.output.dir, "analysis.csv"); colors)
    Makie.save("output/Plummer/Momentum-$(mode).png", fig)

    fig, df = plot_momentum_angular(joinpath(sim.config.output.dir, "analysis.csv"); colors)
    Makie.save("output/Plummer/MomentumAngular-$(mode).png", fig)

    # Plot radii
    FigScale, FigLagrange, df = plot_radii(
        sim.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2();
        colors, times = collect(0.0:0.005:1.0)*u"Gyr",
    )
    Makie.save("output/Plummer/ScaleRadius-$(mode).png", FigScale)
    Makie.save("output/Plummer/LagrangianRadii-$(mode).png", FigLagrange)
    mv("radii.csv", "output/Plummer/radii-$(mode).csv", force = true)
    return true
end

function plot_plummer_unitless(sim::Simulation, mode)
    
end

@testset "Plummer Unitful" begin
    d = generate(PlummerStarCluster(NumSamples = 1000))

    SoftLen = 0.01u"kpc"
    ForceSofteningTable = [SoftLen for i in 1:6]

    TimeEnd = 1.0u"Gyr"
    TimeStep = 0.0003u"Gyr"
    TimeBetweenSnapshots = 0.005u"Gyr"

    analysers = Dict(
        "potential" => total_potential,
        "kinetic" => total_kinetic,
        "momentum" => total_momentum,
        "angularmomentum" => total_angular_momentum,
    )

    @testset "DirectSum CPU" begin
        mode = "DirectSum-serial-const"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSum-MAW-const"
        sim = Simulation(
            deepcopy(d);
            pids = procs(),
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSum-MTW-const"
        sim = Simulation(
            deepcopy(d);
            pids = workers(),
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSum-serial-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSum-MAW-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = procs(),
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSum-MTW-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = workers(),
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)
    end
    
    #TODO potential on GPU
    #=
    @testset "DirectSum GPU" begin
        mode = "DirectSumGPU-const"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            device = GPU(),
            Realtime = false,
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "DirectSumGPU-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            device = GPU(),
            Realtime = false,
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)
    end
    =#

    @testset "Tree" begin
        mode = "Tree-serial-const"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "Tree-MAW-const"
        sim = Simulation(
            deepcopy(d);
            pids = procs(),
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "Tree-MTW-const"
        sim = Simulation(
            deepcopy(d);
            pids = workers(),
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "Tree-serial-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "Tree-MAW-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = procs(),
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "Tree-MTW-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = workers(),
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            GravitySolver = Tree(),
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)
    end

    @testset "FDM" begin
        N = 14

        mode = "FDM-CPU-const"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            GravitySolver = FDM(),
            BoundaryCondition = Vacuum(),
            device = CPU(),
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            EnlargeMesh = 3.0,
            Nx = N,
            Ny = N,
            Nz = N,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "FDM-CPU-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            GravitySolver = FDM(),
            BoundaryCondition = Vacuum(),
            device = CPU(),
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            EnlargeMesh = 3.0,
            Nx = N,
            Ny = N,
            Nz = N,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "FDM-GPU-const"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            GravitySolver = FDM(),
            BoundaryCondition = Vacuum(),
            device = GPU(),
            sparse = false,
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            ForceSofteningTable,
            EnlargeMesh = 3.0,
            Nx = N,
            Ny = N,
            Nz = N,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)

        mode = "FDM-GPU-adapt"
        sim = Simulation(
            deepcopy(d);
            pids = [1],
            GravitySolver = FDM(),
            BoundaryCondition = Vacuum(),
            device = GPU(),
            sparse = false,
            analysers,
            TimeEnd,
            TimeStep = 0.0u"Gyr",
            TimeBetweenSnapshots,
            ForceSofteningTable,
            EnlargeMesh = 3.0,
            Nx = N,
            Ny = N,
            Nz = N,
            OutputDir = "output/Plummer/" * mode,
        );
        run(sim)
        @test plot_plummer(sim, mode)
    end
end

@testset "Plummer Unitless" begin
    d = generate(
        PlummerStarCluster(
            NumSamples = 100,
            VirialRadius = 0.01,
            TotalMass = 1.0e5,
        ),
        nothing,
        constants = Constant(nothing, uAstro)
    )


end