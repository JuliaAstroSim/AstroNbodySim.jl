# Test suite for static tools (no integration tests)
# This file is intended to be included from runtests.jl after:
#     using Test
#     @everywhere using AstroNbodySim
#     @everywhere using AstroNbodySim.AstroSimBase
#     @everywhere using Unitful, UnitfulAstro, PhysicalParticles,
#                     PhysicalMeshes, PhysicalTrees
#     @everywhere astro()

@testset "Static Tools" begin

    # ------------------------------------------------------------------
    # Shared constants / units used throughout the testset
    # ------------------------------------------------------------------
    Gconst  = Constant(uAstro).G
    c_si    = 2.99792458e5u"km/s"            # light speed in km/s
    M_bh    = 1.0e6u"Msun"                   # BH mass
    R_orbit = 10.0u"kpc"                     # orbit radius

    # ------------------------------------------------------------------
    # 1. BlackHole tools
    # ------------------------------------------------------------------
    @testset "BlackHole" begin
        # r_g numeric verification: r_g = 2 * G * M / c^2
        rg_val = r_g(Gconst, M_bh, c_si)
        @test rg_val ≈ 2 * Gconst * M_bh / c_si^2
        @test rg_val > 0u"kpc"

        # Aliases must point to the same function / yield same value
        @test radius_gravity       === r_g
        @test radius_schwarzschild === r_g

        # pseudoNewtonianPotential 4-arg form: pot = -G*M/(R - r_g)
        pot4 = pseudoNewtonianPotential(Gconst, M_bh, c_si, R_orbit)
        @test pot4 ≈ - Gconst * M_bh / (R_orbit - rg_val)

        # pseudoNewtonianAcc 5-arg form with PVector
        n_vec = PVector(1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        acc4  = pseudoNewtonianAcc(Gconst, M_bh, c_si, R_orbit, n_vec)
        @test acc4 isa PVector
        @test acc4 ≈ (- Gconst * M_bh / (R_orbit - rg_val)^2) * n_vec
        @test dot(acc4, n_vec) < 0.0u"kpc^3/Gyr^2"   # unitful compare
    end

    # ------------------------------------------------------------------
    # 2. Elliptic Orbit tools
    # ------------------------------------------------------------------
    @testset "EllipticOrbit" begin
        # Circular orbit: a == r, v = sqrt(G*M/r)
        a_circ = ellipticSemiMajor(Gconst, M_bh, R_orbit,
                                   sqrt(Gconst * M_bh / R_orbit))
        @test a_circ ≈ R_orbit

        # Kepler's third law: T = 2π * sqrt(a^3 / (G*M))
        a_test = 5.0u"kpc"
        T_test = ellipticPeriod(Gconst, M_bh, a_test)
        @test T_test ≈ 2 * pi * sqrt(a_test^3 / (Gconst * M_bh))
        # Period grows as a^{3/2}
        T2 = ellipticPeriod(Gconst, M_bh, 4 * a_test)
        @test T2 / T_test ≈ sqrt(4^3)

        # Scalar eccentricity: e = 1 - v*h / (G*M)
        # For a circular orbit at radius a_test, h = a_test * v_test
        v_test = sqrt(Gconst * M_bh / a_test)
        h_test = a_test * v_test
        e_scalar = eccentricity(Gconst, M_bh, v_test, h_test)
        @test isapprox(e_scalar, 0.0; atol = 1e-10)

        # Vector eccentricity for a circular orbit at radius a_test -> zero vector
        # Must use the same radius a_test as v_test (otherwise it's an elliptical orbit)
        r_vec = PVector(a_test, 0.0u"kpc", 0.0u"kpc")
        v_vec = PVector(0.0u"kpc/Gyr", v_test, 0.0u"kpc/Gyr")
        e_vec = eccentricity(Gconst, M_bh, r_vec, v_vec)
        @test norm(e_vec) < 1.0e-6
    end

    # ------------------------------------------------------------------
    # 3. Time Scales
    # ------------------------------------------------------------------
    @testset "TimeScales" begin
        M_sys = 1.0e10u"Msun"
        R_sys = 10.0u"kpc"
        v_sys = sqrt(Gconst * M_sys / R_sys)
        H0    = 70.0u"km/s/Mpc"
        N_sys = 1000
        rho_sys = 3 * M_sys / (4 * pi * R_sys^3)

        @test typicalvelocity(M_sys, R_sys, Gconst) ≈ sqrt(Gconst * M_sys / R_sys)
        @test crosstime(R_sys, v_sys)               ≈ R_sys / v_sys
        @test hubbletime(H0)                        ≈ 1 / H0
        @test relaxtime(R_sys, v_sys, N_sys)        ≈
              N_sys / (10 * log10(N_sys)) * (R_sys / v_sys)
        @test interactiontime(R_sys, v_sys, N_sys)  ≈ N_sys * (R_sys / v_sys)
        @test dynamicaltime(rho_sys, Gconst)        ≈ sqrt(3 * pi / (16 * Gconst * rho_sys))
        @test freefalltime(rho_sys, Gconst)         ≈ sqrt(3 * pi / (32 * Gconst * rho_sys))
        @test orbitaltime(rho_sys, Gconst)          ≈ sqrt(3 * pi / (Gconst * rho_sys))
    end

    # ------------------------------------------------------------------
    # 4. SofteningLength
    # ------------------------------------------------------------------
    @testset "SofteningLength" begin
        # Construct a 2-particle simulation (Serial on master).
        SoftLenVal = 0.001u"kpc"
        TwoParticleData = StructArray(Star(uAstro, id = i) for i in 1:2)
        TwoParticleData.Pos[1] = PVector(+0.5u"kpc", 0.0u"kpc", 0.0u"kpc")
        TwoParticleData.Pos[2] = PVector(-0.5u"kpc", 0.0u"kpc", 0.0u"kpc")
        TwoParticleData.Mass .= 1.0e10u"Msun"
        ForceTable = [SoftLenVal for _ in 1:6]

        sim = Simulation(deepcopy(TwoParticleData);
                         pids = [1],
                         ForceSofteningTable = ForceTable)
        table = sim.config.grav.ForceSofteningTable  # MVector{6}
        p1 = get_all_data(sim)[1]

        # --- softlen: 4 overloads ---
        @test softlen(STAR,      table) == SoftLenVal          # (c, table)
        @test softlen(p1,        table) == SoftLenVal          # (p, table)
        @test softlen(STAR,      sim)   == SoftLenVal          # (c, sim)
        @test softlen(p1,        sim)   == SoftLenVal          # (p, sim)

        # --- suggest_softlen: 4 overloads ---
        @test suggest_softlen(8.0u"kpc^3", 8) ≈ 1.0u"kpc"                      # (V, N)
        s_data = suggest_softlen(get_all_data(sim))
        @test s_data ≈ suggest_softlen(get_all_data(sim), STAR)                # (data) == (data, c)
        @test suggest_softlen(sim) ≈ s_data                                    # (sim)
        # collection NOT present -> warn + nothing
        # Use 1:2 (positive range) so StructArray can infer the Star type;
        # the actual filter for BLACKHOLE is still empty because the stars are STAR.
        empty_data = StructArray(Star(uAstro, id = i) for i in 1:2)
        s_missing = nothing
        @test_logs (:warn,) s_missing = suggest_softlen(empty_data, BLACKHOLE)
        @test s_missing === nothing

        # --- set_softlen_local!: 2 overloads ---
        new_h_one    = 0.005u"kpc"
        new_h_uniform = 0.002u"kpc"
        # set_softlen_local! is not exported; qualify.
        AstroNbodySim.set_softlen_local!(sim, STAR, new_h_one)                               # (sim, c, h)
        @test sim.config.grav.ForceSofteningTable[Int(STAR)] == new_h_one
        AstroNbodySim.set_softlen_local!(sim, new_h_uniform)                                # (sim, h)
        @test all(==(new_h_uniform), sim.config.grav.ForceSofteningTable)

        # --- set_softlen!: dispatch branches ---
        # Default CPU + DirectSum/Tree -> bcast via set_softlen_local!
        AstroNbodySim.set_softlen_local!(sim, SoftLenVal)
        set_softlen!(sim, new_h_one)
        @test sim.config.grav.ForceSofteningTable[Int(STAR)] == new_h_one
        # CPU + FFT/FDM path -> direct local set
        AstroNbodySim.set_softlen_local!(sim, SoftLenVal)
        AstroNbodySim.set_softlen!(sim, new_h_uniform, FFT(), CPU())
        @test all(==(new_h_uniform), sim.config.grav.ForceSofteningTable)

        # --- suggest_softlen!: updates table to suggested value ---
        AstroNbodySim.set_softlen_local!(sim, SoftLenVal)
        suggested = suggest_softlen(sim)
        suggest_softlen!(sim)
        @test sim.config.grav.ForceSofteningTable[Int(STAR)] ≈ suggested
    end

    # ------------------------------------------------------------------
    # 5. Gadget2 makefile / param IO
    # ------------------------------------------------------------------
    @testset "Gadget2 IO" begin
        tmp = mktempdir()

        # --- write_gadget2_makefile (default + custom options) ---
        mk_path = joinpath(tmp, "Makefile")
        write_gadget2_makefile(mk_path)
        mk_text = read(mk_path, String)
        @test occursin("EXEC   = Gadget2", mk_text)
        @test occursin("OPT += -DPEANOHILBERT", mk_text)
        @test occursin("CC       = mpicc", mk_text)

        write_gadget2_makefile(joinpath(tmp, "Makefile.custom"),
                               "MYOPT1"; compiler = "gcc",
                                          hdf5_lib   = "-L/opt/hdf5/lib")
        @test occursin("CC       = gcc", read(joinpath(tmp, "Makefile.custom"), String))
        @test occursin("OPT += -DMYOPT1", read(joinpath(tmp, "Makefile.custom"), String))

        # --- write_gadget2_param (basic 2-arg form) ---
        g2_path = joinpath(tmp, "gadget2.param")
        write_gadget2_param(g2_path, "ic.dat")
        @test occursin("InitCondFile  	   ic.dat", read(g2_path, String))
        @test occursin("SofteningGas", read(g2_path, String))

        # --- alter_param ---
        alter_param(g2_path, "InitCondFile", "new_ic.dat")
        @test occursin("InitCondFile  new_ic.dat", read(g2_path, String))

        # --- alter_line : default (first match only) vs all = true ---
        alt_path = joinpath(tmp, "alt.param")
        write(alt_path, "Alpha 1\nBeta 2\nAlpha 3\n")
        alter_line(alt_path, "Alpha", "Alpha 100")
        @test occursin("Alpha 100", read(alt_path, String))
        @test occursin("Alpha 3", read(alt_path, String))      # second match preserved

        write(alt_path, "Delta 1\nDelta 2\nDelta 3\n")
        alter_line(alt_path, "Delta", "D 99"; all = true)
        @test occursin("D 99", read(alt_path, String))
        @test !occursin("Delta ", read(alt_path, String))

        # --- add_line : prepends a new line at the top ---
        write(joinpath(tmp, "add.param"), "LineA\nLineB\nLineC\n")
        add_line(joinpath(tmp, "add.param"), "Header")
        @test startswith(read(joinpath(tmp, "add.param"), String), "Header")

        # --- write_gadget4_param 4 overloads (not exported -> qualify) ---
        # (1) keyword form
        g4_path = joinpath(tmp, "gadget4.param")
        AstroNbodySim.write_gadget4_param(g4_path, "ic4.dat";
                                          TimeBegin = 0.0, TimeMax = 0.05)
        @test occursin("InitCondFile", read(g4_path, String))
        @test occursin("SofteningComovingClass0", read(g4_path, String))

        # (2) SimConfig form
        cfg_path = joinpath(tmp, "gadget4_cfg.param")
        FST = [0.01 for _ in 1:6]
        cfg = SimConfig(;
            TimeBegin             = 0.0u"Gyr",
            TimeEnd               = 0.5u"Gyr",
            TimeBetweenSnapshots  = 0.01u"Gyr",
            TimeStep              = 0.0u"Gyr",
            ForceSofteningTable   = FST)
        AstroNbodySim.write_gadget4_param(cfg_path, "ic4_cfg.dat", cfg)
        @test occursin("ErrTolIntAccuracy", read(cfg_path, String))

        # (3) sim via (sim, grav, dev) -- tree solver populates ErrTolTheta
        sim_path  = joinpath(tmp, "gadget4_sim.param")
        sim_path2 = joinpath(tmp, "gadget4_sim2.param")
        TwoParticle = StructArray(Star(uAstro, id = i) for i in 1:2)
        TwoParticle.Pos[1] = PVector(+0.5u"kpc", 0.0u"kpc", 0.0u"kpc")
        TwoParticle.Pos[2] = PVector(-0.5u"kpc", 0.0u"kpc", 0.0u"kpc")
        TwoParticle.Mass .= 1.0e10u"Msun"
        sim_t = Simulation(deepcopy(TwoParticle);
                           pids = [1],
                           GravitySolver = Tree(),
                           ForceSofteningTable = FST,
                           TimeBegin = 0.0u"Gyr",
                           TimeEnd   = 0.2u"Gyr",
                           TimeBetweenSnapshots = 0.01u"Gyr",
                           TimeStep  = 0.0u"Gyr")
        AstroNbodySim.write_gadget4_param(sim_path, "ic4_sim.dat", sim_t,
                                         sim_t.config.solver.grav,
                                         sim_t.config.device.type)
        @test occursin("ErrTolTheta", read(sim_path, String))
        @test occursin("ErrTolForceAcc", read(sim_path, String))

        # (4) sim via the single-argument dispatch -> picks Tree/CPU internally
        AstroNbodySim.write_gadget4_param(sim_path2, "ic4_sim2.dat", sim_t)
        @test occursin("InitCondFile", read(sim_path2, String))
    end

    # ------------------------------------------------------------------
    # 6. SimConfig validation
    # ------------------------------------------------------------------
    @testset "SimConfig Validation" begin
        # ForceSofteningTable length != 6
        @test_throws ArgumentError SimConfig(; ForceSofteningTable = [0.01u"kpc" for _ in 1:3])
        @test_throws ArgumentError SimConfig(; ForceSofteningTable = [0.01u"kpc" for _ in 1:7])
        # ForceSofteningTable contains zero
        FST_zero = [0.01u"kpc" for _ in 1:6]
        FST_zero[3] = 0.0u"kpc"
        @test_throws ArgumentError SimConfig(; ForceSofteningTable = FST_zero)
        # EnlargeMesh negative / zero
        @test_throws ArgumentError SimConfig(; EnlargeMesh = -1.0)
        @test_throws ArgumentError SimConfig(; EnlargeMesh = 0.0)
        # TimeEnd <= TimeBegin
        @test_throws ArgumentError SimConfig(; TimeBegin = 1.0u"Gyr", TimeEnd = 1.0u"Gyr")
        @test_throws ArgumentError SimConfig(; TimeBegin = 2.0u"Gyr", TimeEnd = 1.0u"Gyr")
        # TimeBetweenSnapshots = 0
        @test_throws ArgumentError SimConfig(; TimeBetweenSnapshots = 0.0u"Gyr")
        # SaveRestartFreq <= 0
        @test_throws ArgumentError SimConfig(; SaveRestartFreq = 0.0)
        @test_throws ArgumentError SimConfig(; SaveRestartFreq = -1.5)
        # Sanity: a valid config builds without throwing
        @test_nowarn SimConfig(;
            TimeBegin = 0.0u"Gyr", TimeEnd = 0.1u"Gyr",
            TimeBetweenSnapshots = 0.01u"Gyr", TimeStep = 0.0u"Gyr")
    end

    # ------------------------------------------------------------------
    # 7. tabstring
    # ------------------------------------------------------------------
    @testset "tabstring" begin
        @test AstroNbodySim.tabstring("a\nb")           == "    a\n    b"
        @test AstroNbodySim.tabstring("a\nb", ">>")    == ">>a\n>>b"
        @test AstroNbodySim.tabstring("solo")          == "    solo"
    end

    # ------------------------------------------------------------------
    # 8. LogInfo setproperty! : dict replace semantics
    # ------------------------------------------------------------------
    @testset "LogInfo setproperty!" begin
        li = LogInfo()
        # Assigning the same object is a no-op
        li.analysers = li.analysers
        @test li.analysers === li.analysers
        # Assigning a new Dict: clear + merge
        li.analysers = Dict("FOO" => () -> nothing)
        @test haskey(li.analysers, "FOO")
        @test !haskey(li.analysers, "BAR")
        # timers can also be reassigned
        li3 = LogInfo()
        # setproperty! only mutates the dict in place; `timing` vector retains
        # the length set at construction (8 from DefaultTimer).
        li3.timers = Dict("CUSTOM" => 1)
        @test length(li3.timers) == 1
        @test li3.timing == [UInt64(0) for _ in 1:length(DefaultTimer)]
    end

    # ------------------------------------------------------------------
    # 9. show methods : sprint(show, ...) must not throw
    # ------------------------------------------------------------------
    @testset "show methods" begin
        # OutputConfig is defined in Config.jl but not exported; qualify.
        oc  = AstroNbodySim.OutputConfig()
        gc  = AstroNbodySim.GravityConfig(uAstro)
        sc  = AstroNbodySim.SolverConfig()
        dc  = AstroNbodySim.DeviceConfig()
        tc  = AstroNbodySim.TimeConfig(; units = uAstro,
                                        TimeBegin = 0.0u"Gyr",
                                        TimeEnd = 0.1u"Gyr",
                                        TimeStep = 0.0u"Gyr")
        cfg = SimConfig(; TimeBegin = 0.0u"Gyr", TimeEnd = 0.1u"Gyr",
                          TimeBetweenSnapshots = 0.01u"Gyr",
                          TimeStep = 0.0u"Gyr")
        @test_nowarn sprint(show, oc)
        @test_nowarn sprint(show, gc)
        @test_nowarn sprint(show, sc)
        @test_nowarn sprint(show, dc)
        @test_nowarn sprint(show, tc)
        @test_nowarn sprint(show, cfg)
        @test_nowarn sprint(show, ConstantTimestep(0.001u"Gyr"))
        @test_nowarn sprint(show, LogInfo())
    end

end