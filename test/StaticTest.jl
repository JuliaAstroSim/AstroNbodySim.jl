@testset "Static Test" begin
    # Test data
    
    ## Binary
    G = Constant(uAstro).G
    GUnitless = ustrip(G)

    TimeStep = 0.0001u"Gyr"
    TimeStepUnitless = 0.0001

    ErrTolTimestep = 0.025
    SoftLen = 0.0001u"kpc"
    SoftLenUnitless = 0.0001

    Constants = Constant(uAstro)
    ConstantsUnitless = ustrip(Constants)

    function binary_orbit_period(mass, D, G)
        return sqrt(2pi * D^3 / G / mass)
    end

    function binary_rot_vel(mass, D, G)
        return sqrt(G * mass / D / 2.0)
    end

    function binary_potential(mass, D, G, h)
        # 0.5 * \sum(Pot_i)
        return -1.0 * G * mass^2 / sqrt(D^2 + h^2) #+ -1.0 * G * mass^2 / sqrt(h^2)
    end

    function binary_kinetic(mass, D, G)
        v = binary_rot_vel(mass, D, G)
        return mass * v^2
    end

    function binary_energy(mass, D, G, h)
        p = binary_potential(mass, D, G, h)
        k = binary_kinetic(mass, D, G)
        return p + k
    end

    function binary_momentum(mass, vel)
        return 0.0 * mass * vel * PVector()
    end

    function binary_momentum_angular(mass, vel, D)
        return mass * vel * D * PVector(0.0, 0.0, 1.0)
    end

    function binary_force(mass, D, G, h)
        return G * mass / (D^2 + h^2) / sqrt((D^2 + h^2)) * D
    end

    BinaryMass = 1.0e10u"Msun"
    BinaryD = 2.0u"kpc"

    BinaryVel = binary_rot_vel(BinaryMass, BinaryD, G)
    BinaryAcc = binary_force(BinaryMass, BinaryD, G, SoftLen)
    BinaryPotential = binary_potential(BinaryMass, BinaryD, G, SoftLen)
    BinaryKinetic = binary_kinetic(BinaryMass, BinaryD, G)
    BinaryEnergy = binary_energy(BinaryMass, BinaryD, G, SoftLen)
    BinaryMomentum = binary_momentum(BinaryMass, BinaryVel)
    BinaryMomentumAngular = binary_momentum_angular(BinaryMass, BinaryVel, BinaryD)

    BinaryData = StructArray(Star(uAstro, id = i) for i in 1:2)
    BinaryData.Pos[1] = PVector(+0.5BinaryD, 0.0u"kpc", 0.0u"kpc")
    BinaryData.Pos[2] = PVector(-0.5BinaryD, 0.0u"kpc", 0.0u"kpc")
    BinaryData.Vel[1] = PVector(0.0u"kpc/Gyr", +BinaryVel, 0.0u"kpc/Gyr")
    BinaryData.Vel[2] = PVector(0.0u"kpc/Gyr", -BinaryVel, 0.0u"kpc/Gyr")
    BinaryData.Mass .= BinaryMass



    BinaryMassUnitless = 1.0e10
    BinaryDUnitless = 2.0

    BinaryVelUnitless = binary_rot_vel(BinaryMassUnitless, BinaryDUnitless, GUnitless)
    BinaryAccUnitless = binary_force(BinaryMassUnitless, BinaryDUnitless, GUnitless, SoftLenUnitless)
    BinaryPotentialUnitless = binary_potential(BinaryMassUnitless, BinaryDUnitless, GUnitless, SoftLenUnitless)
    BinaryKineticUnitless = binary_kinetic(BinaryMassUnitless, BinaryDUnitless, GUnitless)
    BinaryEnergyUnitless = binary_energy(BinaryMassUnitless, BinaryDUnitless, GUnitless, SoftLenUnitless)
    BinaryMomentumUnitless = binary_momentum(BinaryMassUnitless, BinaryVelUnitless)
    BinaryMomentumAngularUnitless = binary_momentum_angular(BinaryMassUnitless, BinaryVelUnitless, BinaryDUnitless)

    BinaryDataUnitless = StructArray(Star(id = i) for i in 1:2);
    BinaryDataUnitless.Pos[1] = PVector(+0.5BinaryDUnitless, 0.0, 0.0)
    BinaryDataUnitless.Pos[2] = PVector(-0.5BinaryDUnitless, 0.0, 0.0)
    BinaryDataUnitless.Vel[1] = PVector(0.0, +BinaryVelUnitless, 0.0)
    BinaryDataUnitless.Vel[2] = PVector(0.0, -BinaryVelUnitless, 0.0)
    BinaryDataUnitless.Mass .= BinaryMassUnitless

    ## Box
    BoxPos = [
        PVector(+0.1 + 1.9, +0.1 + 1.9, +0.1 + 1.9, u"kpc"),
        PVector(-0.1 + 1.9, +0.1 + 1.9, +0.1 + 1.9, u"kpc"),
        PVector(+0.1 + 1.9, -0.1 + 1.9, +0.1 + 1.9, u"kpc"),
        PVector(+0.1 + 1.9, +0.1 + 1.9, -0.1 + 1.9, u"kpc"),
        PVector(-0.1 + 1.9, -0.1 + 1.9, +0.1 + 1.9, u"kpc"),
        PVector(-0.1 + 1.9, +0.1 + 1.9, -0.1 + 1.9, u"kpc"),
        PVector(+0.1 + 1.9, -0.1 + 1.9, -0.1 + 1.9, u"kpc"),
        PVector(-0.1 + 1.9, -0.1 + 1.9, -0.1 + 1.9, u"kpc"),
        PVector(+1.0, +1.0, +1.0, u"kpc"),
        PVector(-1.0, -1.0, -1.0, u"kpc"),
    ]
    BoxData = StructArray(Star(uAstro, id = i) for i in 1:10);
    assign_particles(BoxData, :Pos, BoxPos)
    assign_particles(BoxData, :Mass, 1.0e9u"Msun")

    BoxPosUnitless = [
        PVector(+0.1 + 1.9, +0.1 + 1.9, +0.1 + 1.9),
        PVector(-0.1 + 1.9, +0.1 + 1.9, +0.1 + 1.9),
        PVector(+0.1 + 1.9, -0.1 + 1.9, +0.1 + 1.9),
        PVector(+0.1 + 1.9, +0.1 + 1.9, -0.1 + 1.9),
        PVector(-0.1 + 1.9, -0.1 + 1.9, +0.1 + 1.9),
        PVector(-0.1 + 1.9, +0.1 + 1.9, -0.1 + 1.9),
        PVector(+0.1 + 1.9, -0.1 + 1.9, -0.1 + 1.9),
        PVector(-0.1 + 1.9, -0.1 + 1.9, -0.1 + 1.9),
        PVector(+1.0, +1.0, +1.0),
        PVector(-1.0, -1.0, -1.0),
    ]
    BoxDataUnitless = StructArray(Star(id = i) for i in 1:10);
    assign_particles(BoxDataUnitless, :Pos, BoxPosUnitless)
    assign_particles(BoxDataUnitless, :Mass, 1.0e9)

    zerovalues = ZeroValue(uAstro)
    zerovaluesUnitless = ZeroValue(nothing)

    @testset "DirectSum" begin
        @testset "Binary" begin
            BinarySim_SerialOnMaster = Simulation(
                deepcopy(BinaryData),
                pids = [1];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_SerialOnWorker = Simulation(
                deepcopy(BinaryData),
                pids = [2];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_MasterAsWorker = Simulation(
                deepcopy(BinaryData),
                pids = [1,3];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_MasterToWorker = Simulation(
                deepcopy(BinaryData),
                pids = [2,4];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySims = [
                BinarySim_SerialOnMaster,
                BinarySim_SerialOnWorker,
                BinarySim_MasterAsWorker,
                BinarySim_MasterToWorker,
            ];
    
            compute_force.(BinarySims)
    
            @testset "Force" begin
                for sim in BinarySims
                    data = get_all_data(sim)
                    a1 = data.Acc[1].x
                    a2 = data.Acc[2].x
                    a0 = a1 + a2
                    @test a0 ≈ zero(a0)
                    @test abs(a1) ≈ BinaryAcc
                    @test abs(a2) ≈ BinaryAcc
                end
            end
    
            compute_potential.(BinarySims)
    
            @testset "Energy" begin
                for sim in BinarySims
                    p = total_potential(sim)
                    k = total_kinetic(sim)
                    e = total_energy(sim)
                    @test p ≈ BinaryPotential
                    @test k ≈ BinaryKinetic
                    @test e ≈ BinaryEnergy
                end
            end
    
            AstroNbodySim.init_timesteps.(BinarySims)
    
            @testset "Timestep" begin
                for sim in BinarySims
                    p = first(get_all_data(sim))
                    h = sim.config.grav.ForceSofteningTable[Int(p.Collection)]
                    dt_theoretical = sqrt(2.0 * sim.config.time.step.ErrTolTimestep * h / norm(p.Acc))
                    dt_sim = (p.Ti_endstep - p.Ti_begstep) * sim.config.time.step.TimeInterval
                    @test dt_sim < dt_theoretical
                end
            end
    
            @testset "Momentum" begin
                for sim in BinarySims
                    P = total_momentum(sim)
                    L = total_angular_momentum(sim)
                    @test P ≈ BinaryMomentum
                    @test L ≈ BinaryMomentumAngular
                end
            end
        end
    
        @testset "Binary Unitless" begin
            BinarySimUnitless_SerialOnMaster = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [1];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_SerialOnWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [2];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_MasterAsWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [1,3];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_MasterToWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [2,4];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimsUnitless = [
                BinarySimUnitless_SerialOnMaster,
                BinarySimUnitless_SerialOnWorker,
                BinarySimUnitless_MasterAsWorker,
                BinarySimUnitless_MasterToWorker,
            ];
    
            compute_force.(BinarySimsUnitless)
    
            @testset "Force" begin
                for sim in BinarySimsUnitless
                    data = get_all_data(sim)
                    a1 = data.Acc[1].x
                    a2 = data.Acc[2].x
                    a0 = a1 + a2
                    @test a0 ≈ zero(a0)
                    @test abs(a1) ≈ BinaryAccUnitless
                    @test abs(a2) ≈ BinaryAccUnitless
                end
            end
    
            compute_potential.(BinarySimsUnitless)
    
            @testset "Energy" begin
                for sim in BinarySimsUnitless
                    p = total_potential(sim)
                    k = total_kinetic(sim)
                    e = total_energy(sim)
                    @test p ≈ BinaryPotentialUnitless
                    @test k ≈ BinaryKineticUnitless
                    @test e ≈ BinaryEnergyUnitless
                end
            end
    
            AstroNbodySim.init_timesteps.(BinarySimsUnitless)
    
            @testset "Timestep" begin
                for sim in BinarySimsUnitless
                    p = first(get_all_data(sim))
                    h = sim.config.grav.ForceSofteningTable[Int(p.Collection)]
                    dt_theoretical = sqrt(2.0 * sim.config.time.step.ErrTolTimestep * h / norm(p.Acc))
                    dt_sim = (p.Ti_endstep - p.Ti_begstep) * sim.config.time.step.TimeInterval
                    @test dt_sim < dt_theoretical
                end
            end
    
            @testset "Momentum" begin
                for sim in BinarySimsUnitless
                    P = total_momentum(sim)
                    L = total_angular_momentum(sim)
                    @test P ≈ BinaryMomentumUnitless
                    @test L ≈ BinaryMomentumAngularUnitless
                end
            end
        end
    
        @testset "Box" begin
            BoxSim_SerialOnMaster = Simulation(
                deepcopy(BoxData),
                pids = [1];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_SerialOnWorker = Simulation(
                deepcopy(BoxData),
                pids = [2];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_MasterAsWorker = Simulation(
                deepcopy(BoxData),
                pids = [1,3];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_MasterToWorker = Simulation(
                deepcopy(BoxData),
                pids = [2,4];
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSims = [
                BoxSim_SerialOnMaster,
                BoxSim_SerialOnWorker,
                BoxSim_MasterAsWorker,
                BoxSim_MasterToWorker,
            ];
    
            compute_force.(BoxSims)
    
            @testset "Force" begin
                pos = PVector(u"kpc")
                BoxAcc = AstroNbodySim.compute_force_at_point(pos, BoxData, G, SoftLen)
                for sim in BoxSims
                    @test first(compute_force(sim, pos, SoftLen)) ≈ BoxAcc
                end
            end
    
            compute_potential.(BoxSims)
    
            @testset "Energy" begin
                pos = PVector(u"kpc")
                BoxPot = AstroNbodySim.compute_unit_potential_at_point(pos, BoxData, G, SoftLen)
                for sim in BoxSims
                    @test first(compute_potential(sim, [pos], SoftLen, sim.config.solver.grav, sim.config.device.type)) ≈ BoxPot
                end
            end
        end
    
        @testset "Box Unitless" begin
            BoxSimUnitless_SerialOnMaster = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [1];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_SerialOnWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [2];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_MasterAsWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [1,3];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_MasterToWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [2,4];
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimsUnitless = [
                BoxSimUnitless_SerialOnMaster,
                BoxSimUnitless_SerialOnWorker,
                BoxSimUnitless_MasterAsWorker,
                BoxSimUnitless_MasterToWorker,
            ];
    
            compute_force.(BoxSimsUnitless)
    
            @testset "Force" begin
                pos = PVector()
                BoxAccUnitless = AstroNbodySim.compute_force_at_point(pos, BoxDataUnitless, GUnitless, SoftLenUnitless)
                for sim in BoxSimsUnitless
                    @test first(compute_force(sim, pos, SoftLenUnitless)) ≈ BoxAccUnitless
                end
            end
    
            compute_potential.(BoxSimsUnitless)
    
            @testset "Energy" begin
                pos = PVector()
                BoxPotUnitless = AstroNbodySim.compute_unit_potential_at_point(pos, BoxDataUnitless, GUnitless, SoftLenUnitless)
                for sim in BoxSimsUnitless
                    @test first(compute_potential(sim, [pos], SoftLenUnitless, sim.config.solver.grav, sim.config.device.type)) ≈ BoxPotUnitless
                end
            end
        end
    end

    @testset "Tree" begin
        @testset "Binary" begin
            BinarySim_SerialOnMaster = Simulation(
                deepcopy(BinaryData),
                pids = [1];
                ErrTolTimestep,
                GravitySolver = Tree(),
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_SerialOnWorker = Simulation(
                deepcopy(BinaryData),
                pids = [2];
                ErrTolTimestep,
                GravitySolver = Tree(),
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_MasterAsWorker = Simulation(
                deepcopy(BinaryData),
                pids = [1,3];
                ErrTolTimestep,
                GravitySolver = Tree(),
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySim_MasterToWorker = Simulation(
                deepcopy(BinaryData),
                pids = [2,4];
                ErrTolTimestep,
                GravitySolver = Tree(),
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BinarySims = [
                BinarySim_SerialOnMaster,
                BinarySim_SerialOnWorker,
                BinarySim_MasterAsWorker,
                BinarySim_MasterToWorker,
            ];
    
            compute_force.(BinarySims)
    
            @testset "Force" begin
                for sim in BinarySims
                    data = get_all_data(sim)
                    a1 = data.Acc[1].x
                    a2 = data.Acc[2].x
                    a0 = a1 + a2
                    @test a0 ≈ zero(a0)
                    @test abs(a1) ≈ BinaryAcc
                    @test abs(a2) ≈ BinaryAcc
                end
            end
    
            compute_potential.(BinarySims)
    
            @testset "Energy" begin
                for sim in BinarySims
                    p = total_potential(sim)
                    k = total_kinetic(sim)
                    e = total_energy(sim)
                    @test p ≈ BinaryPotential
                    @test k ≈ BinaryKinetic
                    @test e ≈ BinaryEnergy
                end
            end
    
            AstroNbodySim.init_timesteps.(BinarySims)
    
            @testset "Timestep" begin
                for sim in BinarySims
                    p = first(get_all_data(sim))
                    h = sim.config.grav.ForceSofteningTable[Int(p.Collection)]
                    dt_theoretical = sqrt(2.0 * sim.config.time.step.ErrTolTimestep * h / norm(p.Acc))
                    dt_sim = (p.Ti_endstep - p.Ti_begstep) * sim.config.time.step.TimeInterval
                    @test dt_sim < dt_theoretical
                end
            end
    
            @testset "Momentum" begin
                for sim in BinarySims
                    P = total_momentum(sim)
                    L = total_angular_momentum(sim)
                    @test P ≈ BinaryMomentum
                    @test L ≈ BinaryMomentumAngular
                end
            end
        end
    
        @testset "Binary Unitless" begin
            BinarySimUnitless_SerialOnMaster = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [1];
                ErrTolTimestep,
                units = nothing,
                GravitySolver = Tree(),
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_SerialOnWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [2];
                ErrTolTimestep,
                units = nothing,
                GravitySolver = Tree(),
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_MasterAsWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [1,3];
                ErrTolTimestep,
                units = nothing,
                GravitySolver = Tree(),
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimUnitless_MasterToWorker = Simulation(
                deepcopy(BinaryDataUnitless),
                pids = [2,4];
                ErrTolTimestep,
                units = nothing,
                GravitySolver = Tree(),
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BinarySimsUnitless = [
                BinarySimUnitless_SerialOnMaster,
                BinarySimUnitless_SerialOnWorker,
                BinarySimUnitless_MasterAsWorker,
                BinarySimUnitless_MasterToWorker,
            ];
    
            compute_force.(BinarySimsUnitless)
    
            @testset "Force" begin
                for sim in BinarySimsUnitless
                    data = get_all_data(sim)
                    a1 = data.Acc[1].x
                    a2 = data.Acc[2].x
                    a0 = a1 + a2
                    @test a0 ≈ zero(a0)
                    @test abs(a1) ≈ BinaryAccUnitless
                    @test abs(a2) ≈ BinaryAccUnitless
                end
            end
    
            compute_potential.(BinarySimsUnitless)
    
            @testset "Energy" begin
                for sim in BinarySimsUnitless
                    p = total_potential(sim)
                    k = total_kinetic(sim)
                    e = total_energy(sim)
                    @test p ≈ BinaryPotentialUnitless
                    @test k ≈ BinaryKineticUnitless
                    @test e ≈ BinaryEnergyUnitless
                end
            end
    
            AstroNbodySim.init_timesteps.(BinarySimsUnitless)
    
            @testset "Timestep" begin
                for sim in BinarySimsUnitless
                    p = first(get_all_data(sim))
                    h = sim.config.grav.ForceSofteningTable[Int(p.Collection)]
                    dt_theoretical = sqrt(2.0 * sim.config.time.step.ErrTolTimestep * h / norm(p.Acc))
                    dt_sim = (p.Ti_endstep - p.Ti_begstep) * sim.config.time.step.TimeInterval
                    @test dt_sim < dt_theoretical
                end
            end
    
            @testset "Momentum" begin
                for sim in BinarySimsUnitless
                    P = total_momentum(sim)
                    L = total_angular_momentum(sim)
                    @test P ≈ BinaryMomentumUnitless
                    @test L ≈ BinaryMomentumAngularUnitless
                end
            end
        end
    
        @testset "Box" begin
            BoxSim_SerialOnMaster = Simulation(
                deepcopy(BoxData),
                pids = [1];
                ErrTolTimestep,
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_SerialOnWorker = Simulation(
                deepcopy(BoxData),
                pids = [2];
                ErrTolTimestep,
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_MasterAsWorker = Simulation(
                deepcopy(BoxData),
                pids = [1,3];
                ErrTolTimestep,
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSim_MasterToWorker = Simulation(
                deepcopy(BoxData),
                pids = [2,4];
                ErrTolTimestep,
                ForceSofteningTable = [SoftLen for i in 1:6],
            );
    
            BoxSims = [
                BoxSim_SerialOnMaster,
                BoxSim_SerialOnWorker,
                BoxSim_MasterAsWorker, #TODO
                BoxSim_MasterToWorker,
            ];
    
            compute_force.(BoxSims)
    
            @testset "Force" begin
                pos = PVector(u"kpc")
                BoxAcc = AstroNbodySim.compute_force_at_point(pos, BoxData, G, SoftLen)
                for sim in BoxSims
                    @test first(compute_force(sim, pos, SoftLen)) ≈ BoxAcc
                end
            end
    
            compute_potential.(BoxSims)
    
            @testset "Energy" begin
                pos = PVector(u"kpc")
                BoxPot = AstroNbodySim.compute_unit_potential_at_point(pos, BoxData, G, SoftLen)
                for sim in BoxSims
                    @test first(compute_potential(sim, [pos], SoftLen, sim.config.solver.grav, sim.config.device.type)) ≈ BoxPot
                end
            end
        end
    
        @testset "Box Unitless" begin
            BoxSimUnitless_SerialOnMaster = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [1];
                ErrTolTimestep,
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_SerialOnWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [2];
                ErrTolTimestep,
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_MasterAsWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [1,3];
                ErrTolTimestep,
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimUnitless_MasterToWorker = Simulation(
                deepcopy(BoxDataUnitless),
                pids = [2,4];
                ErrTolTimestep,
                units = nothing,
                constants = ConstantsUnitless,
                ForceSofteningTable = [SoftLenUnitless for i in 1:6],
            );
    
            BoxSimsUnitless = [
                BoxSimUnitless_SerialOnMaster,
                BoxSimUnitless_SerialOnWorker,
                BoxSimUnitless_MasterAsWorker,
                BoxSimUnitless_MasterToWorker,
            ];
    
            compute_force.(BoxSimsUnitless)
    
            @testset "Force" begin
                pos = PVector()
                BoxAccUnitless = AstroNbodySim.compute_force_at_point(pos, BoxDataUnitless, GUnitless, SoftLenUnitless)
                for sim in BoxSimsUnitless
                    @test first(compute_force(sim, pos, SoftLenUnitless)) ≈ BoxAccUnitless
                end
            end
    
            compute_potential.(BoxSimsUnitless)
    
            @testset "Energy" begin
                pos = PVector()
                BoxPotUnitless = AstroNbodySim.compute_unit_potential_at_point(pos, BoxDataUnitless, GUnitless, SoftLenUnitless)
                for sim in BoxSimsUnitless
                    @test first(compute_potential(sim, [pos], SoftLenUnitless, sim.config.solver.grav, sim.config.device.type)) ≈ BoxPotUnitless
                end
            end
        end
    end
end