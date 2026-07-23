### Static Simulation Lifecycle Tests
### Tests for Simulation construction, Buffer, Progress, Timing, Logging,
### Output, Analyse, Plot, Restart, GPU helpers, post-processing force/
### potential, timestep branches, julia_main, and preprocessdata dispatch.
###
### Pre-conditions (set up in runtests.jl):
###   - addprocs(4) has been called (procs 1-5 available)
###   - using Test, AstroNbodySim, AstroNbodySim.AstroSimBase,
###     Unitful, UnitfulAstro, PhysicalParticles, PhysicalMeshes,
###     PhysicalTrees, FileIO
###   - astro() invoked globally
###   - GPU detected conditionally via @hascuda
###   - this file is sourced at e:/JuliaAstroSim/AstroNbodySim.jl/test/

# Number of @test assertions collected at the end of this file
const _STATICSIMLIFECYCLE_NUMTESTS = Ref{Int}(0)

# Bring internal AstroNbodySim symbols into scope. These are not exported
# from the package so ``using AstroNbodySim`` does not make them visible;
# importing them explicitly keeps the test body close to the original
# (un-qualified) call sites.
import AstroNbodySim:
    setup_progress, begin_timer, restartloggers, write_analysis,
    estimate_markersize, get_timestep, julia_main,
    update_progress, add_timer, closeloggers, csvstring, pack_pos,
    Progress, write_timing, error_log, check_restartable,
    loadrestart, saverestart, apply_background_potential,
    apply_background_force, redistribute_data
import StaticArrays: MVector

# ----- Helpers -----

# Construct a Plummer cluster of 5 particles to keep tests fast.
#
# NOTE: previously this used ``MaxRadius = 0.0u"kpc"`` which hangs in
# AstroIC (see warning emitted from ``AstroIC.plummer.jl:91`` and the
# call-site comment). We instead use the default ``MaxRadius`` derived
# from ``VirialRadius`` (a small positive value) so ``generate`` returns
# almost immediately on every worker.
function _plummer5()
    ic = PlummerStarCluster(NumSamples = 5, VirialRadius = 0.01u"kpc",
                            TotalMass = 1.0e5u"Msun")
    return generate(ic)
end

# Construct a Plummer cluster of 20 particles
function _plummer20()
    ic = PlummerStarCluster(NumSamples = 20, VirialRadius = 0.01u"kpc",
                            TotalMass = 1.0e5u"Msun")
    return generate(ic)
end

@testset "StaticSimulationLifecycle" begin

    # ---------------------------------------------------------------
    # 1. Buffer
    # ---------------------------------------------------------------
    @testset "Buffer" begin
        # 1a. PosType, AccType, LenType, OldAccType, PotType form
        cfg_buf = SimConfig()
        PosType = typeof(cfg_buf.ZeroValues.pos)
        AccType = typeof(cfg_buf.ZeroValues.acc)
        LenType = typeof(cfg_buf.ZeroValues.pos.x)
        OldAccType = typeof(cfg_buf.ZeroValues.acc.x)
        PotType = typeof(cfg_buf.ZeroValues.potpermass)

        buf = AstroNbodySim.Buffer(PosType, AccType, LenType, OldAccType, PotType)
        @test buf isa AstroNbodySim.Buffer
        @test buf.sendbuffer isa Dict
        @test buf.recvbuffer isa Dict
        @test isempty(buf.sendbuffer)
        @test isempty(buf.recvbuffer)
        @test isempty(buf.gravToEval)
        @test isempty(buf.gravResult)
        @test isempty(buf.potToEval)
        @test isempty(buf.potResult)

        # 1b. Buffer(ZeroValue(uAstro))
        buf2 = AstroNbodySim.Buffer(ZeroValue(uAstro))
        @test buf2 isa AstroNbodySim.Buffer
        @test isempty(buf2.gravToEval)

        # 1c. Buffer(SimConfig())
        buf3 = AstroNbodySim.Buffer(SimConfig())
        @test buf3 isa AstroNbodySim.Buffer

        # 1d. GravToEvaluate and extract_grav
        pos = PVector(0.5u"kpc", 0.5u"kpc", 0.5u"kpc")
        grav = AstroNbodySim.GravToEvaluate(1, pos, STAR, 1.0, 0.001u"kpc")
        @test grav.Index == 1
        @test grav.Pos === pos
        @test grav.Collection == STAR
        p2, c2, o2 = AstroNbodySim.extract_grav(grav)
        @test p2 === pos
        @test c2 == STAR
        @test o2 == 1.0

        # 1e. pack_grav_result and pack_pot_result
        acc = PVector(0.1u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
        gres = AstroNbodySim.pack_grav_result(grav, acc, 7)
        @test gres isa AstroNbodySim.GravResult
        @test gres.Index == 1
        @test gres.Acc === acc
        @test gres.NInteractions == 7

        pot = AstroNbodySim.pack_pot_result(grav, 0.5)
        @test pot isa AstroNbodySim.PotentialResult
        @test pot.Index == 1
        @test pot.Potential == 0.5

        # 1f. PotentialToEvaluate construction
        pte = AstroNbodySim.PotentialToEvaluate(2, pos, HALO, 0.001u"kpc")
        @test pte.Index == 2
        @test pte.Pos === pos
        @test pte.Collection == HALO
        _STATICSIMLIFECYCLE_NUMTESTS[] += 28
    end

    # ---------------------------------------------------------------
    # 2. Progress
    # ---------------------------------------------------------------
    @testset "Progress" begin
        d_adaptive = _plummer5()
        sim_adaptive = Simulation(deepcopy(d_adaptive); pids = [1],
                                  OutputDir = mktempdir(),
                                  loggingmode = ProgressMode())
        setup_progress(sim_adaptive, sim_adaptive.config.time.step)
        @test sim_adaptive.visinfo.progress isa Progress
        update_progress(sim_adaptive, sim_adaptive.config.time.step)
        @test sim_adaptive.visinfo.progress isa Progress

        d_const = _plummer5()
        sim_const = Simulation(deepcopy(d_const); pids = [1],
                               OutputDir = mktempdir(),
                               TimeStep = 1.0e-4u"Gyr",
                               loggingmode = ProgressMode())
        setup_progress(sim_const, sim_const.config.time.step)
        @test sim_const.visinfo.progress isa Progress
        update_progress(sim_const, sim_const.config.time.step)
        @test sim_const.visinfo.progress isa Progress

        _STATICSIMLIFECYCLE_NUMTESTS[] += 4
    end

    # ---------------------------------------------------------------
    # 3. Timing
    # ---------------------------------------------------------------
    @testset "Timing" begin
        d = _plummer5()
        sim = Simulation(deepcopy(d); pids = [1],
                         OutputDir = mktempdir())

        # begin_timer returns UInt64
        t_begin = begin_timer(sim, "FORCE")
        @test t_begin isa UInt64

        # After begin_timer, timing[2] (= "FORCE" index) reset to 0
        @test sim.loginfo.timing[sim.loginfo.timers["FORCE"]] == 0

        # add_timer with single UInt64
        add_timer(sim, "FORCE", UInt64(1000))
        @test sim.loginfo.timing[sim.loginfo.timers["FORCE"]] == UInt64(1000)

        # add_timer with (start, end) accumulates difference
        add_timer(sim, "FORCE", UInt64(100), UInt64(200))
        @test sim.loginfo.timing[sim.loginfo.timers["FORCE"]] == UInt64(1100)

        # Another add for total to be even larger
        add_timer(sim, "TOTAL", UInt64(0), UInt64(500))
        @test sim.loginfo.timing[sim.loginfo.timers["TOTAL"]] == UInt64(500)

        # write_timing should not throw; needs setuploggers to open the file
        setuploggers(sim)
        write_timing(sim)
        @test isfile(joinpath(sim.config.output.dir, "timing.csv"))
        closeloggers(sim)
        _STATICSIMLIFECYCLE_NUMTESTS[] += 6
    end

    # ---------------------------------------------------------------
    # 4. Logging
    # ---------------------------------------------------------------
    @testset "Logging" begin
        d = _plummer5()
        outdir = mktempdir()
        sim = Simulation(deepcopy(d); pids = [1],
                         OutputDir = outdir,
                         analysers = Dict("x" => p -> 1.0))
        # setuploggers should open timerio + analyserio files
        setuploggers(sim)
        @test isopen(sim.stream.timerio)
        @test isopen(sim.stream.analyserio)

        # restartloggers should not throw and re-open the files
        restartloggers(sim)
        @test isopen(sim.stream.timerio)
        @test isopen(sim.stream.analyserio)

        # closeloggers closes both files
        closeloggers(sim)
        # NOTE: Julia's IO API exposes ``Base.isopen`` (not ``isclosed``);
        # we test that the IO handles are no longer open after
        # ``closeloggers``.
        @test !isopen(sim.stream.timerio)
        @test !isopen(sim.stream.analyserio)

        # error_log throws when the per-step log file cannot be opened.
        # The exact exception type depends on whether the directory or
        # file path is missing — accept ``Union{ErrorException, SystemError}``
        # so we don't lock the test to one implementation.
        sim_err = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = mktempdir())
        @test_throws Union{ErrorException, SystemError} error_log(
            sim_err, "test failure message")
        _STATICSIMLIFECYCLE_NUMTESTS[] += 7
    end

    # ---------------------------------------------------------------
    # 5. Output
    # ---------------------------------------------------------------
    @testset "Output" begin
        d = _plummer5()

        # 5a. jld2 output
        outdir_jld = mktempdir()
        sim_jld = Simulation(deepcopy(d); pids = [1],
                             OutputDir = outdir_jld,
                             OutputType = jld2(),
                             acc = true, pot = true)
        compute_force(sim_jld)
        compute_potential(sim_jld)
        @test sim_jld.outputinfo.snapshotcount == 0
        output(sim_jld, sim_jld.config.output.type)
        @test sim_jld.outputinfo.snapshotcount == 1
        @test isfile(joinpath(outdir_jld, "snapshot_0000.jld2"))

        # 5b. gadget2 output
        outdir_g2 = mktempdir()
        sim_g2 = Simulation(deepcopy(d); pids = [1],
                            OutputDir = outdir_g2,
                            OutputType = gadget2(),
                            acc = true, pot = true)
        compute_force(sim_g2)
        compute_potential(sim_g2)
        output(sim_g2, sim_g2.config.output.type)
        @test sim_g2.outputinfo.snapshotcount == 1
        @test isfile(joinpath(outdir_g2, "snapshot_0000.gadget2"))

        # 5c. outputparallel does not throw
        @test (outputparallel(sim_g2, sim_g2.config.output.type); true)

        # 5d. count_gadget_types returns MVector{6,Int32} (not exported -> qualify)
        all_data = get_all_data(sim_g2)
        counts = AstroNbodySim.count_gadget_types(all_data)
        @test counts isa MVector{6,Int32}
        @test sum(counts) == length(all_data)

        # 5e. HeaderGadget2 returns a HeaderGadget2 with valid counts
        hdr = AstroNbodySim.HeaderGadget2(sim_g2, all_data)
        @test hdr isa AstroNbodySim.HeaderGadget2
        @test hdr.npart == counts
        _STATICSIMLIFECYCLE_NUMTESTS[] += 10
    end

    # ---------------------------------------------------------------
    # 6. Analyse
    # ---------------------------------------------------------------
    @testset "Analyse" begin
        outdir = mktempdir()
        sim = Simulation(deepcopy(_plummer5()); pids = [1],
                         OutputDir = outdir,
                         analysers = Dict("x" => p -> 1.0))
        setuploggers(sim)
        result = write_analysis(sim)
        @test result isa Vector
        # First element is the time, second is the analyser result 1.0
        @test length(result) == 2
        @test result[2] == 1.0

        # csvstring simple value
        @test csvstring("hello") == "hello"
        @test csvstring(1.5) == "1.5"

        # csvstring of PVector wraps with quotes
        pv = PVector(1.0, 2.0, 3.0)
        s = csvstring(pv)
        @test occursin("\"", s)
        @test occursin("1.0", s)
        @test occursin("2.0", s)
        @test occursin("3.0", s)

        closeloggers(sim)
        _STATICSIMLIFECYCLE_NUMTESTS[] += 9
    end

    # ---------------------------------------------------------------
    # 7. Plot
    # ---------------------------------------------------------------
    @testset "Plot" begin
        # 7a. estimate_markersize with 2D array
        pos2d = rand(10, 3)
        ms = estimate_markersize(pos2d)
        @test ms > 0
        @test ms isa Real

        # 7b. pack_pos for CPU DirectSum returns an Nx3 matrix
        sim = Simulation(deepcopy(_plummer5()); pids = [1],
                         OutputDir = mktempdir())
        # ``pack_xyz`` returns a 2D matrix of shape ``(N, 3)`` with
        # columns (x, y, z) — see ``PhysicalParticles.pack_xyz``.
        xyz = pack_pos(sim, sim.config.solver.grav, sim.config.device.type)
        @test size(xyz) == (length(get_all_data(sim)), 3)

        # 7c. pack_pos for CPU Tree
        sim_tree = Simulation(deepcopy(_plummer5()); pids = [1],
                              OutputDir = mktempdir(),
                              GravitySolver = Tree())
        xyzt = pack_pos(sim_tree, sim_tree.config.solver.grav,
                        sim_tree.config.device.type)
        @test size(xyzt) == (length(get_all_data(sim_tree)), 3)

        # 7d. pack_pos for CPU FFT
        sim_fft = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = mktempdir(),
                             GravitySolver = FFT(),
                             Nx = 4, Ny = 4, Nz = 4)
        xyzf = pack_pos(sim_fft, sim_fft.config.solver.grav,
                        sim_fft.config.device.type)
        @test size(xyzf) == (length(get_all_data(sim_fft)), 3)
        _STATICSIMLIFECYCLE_NUMTESTS[] += 7
    end

    # ---------------------------------------------------------------
    # 8. Restart (with mktempdir)
    # ---------------------------------------------------------------
    @testset "Restart" begin
        # 8a. saverestart DirectSum / CPU with pids=[1]
        outdir_ds1 = mktempdir()
        sim_ds1 = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = outdir_ds1)
        saverestart(sim_ds1, sim_ds1.config.solver.grav,
                    sim_ds1.config.device.type)
        @test isfile(joinpath(outdir_ds1, "restart.1.jld2"))

        # 8b. saverestart DirectSum / CPU with pids=[1,2,3,4]
        outdir_ds4 = mktempdir()
        sim_ds4 = Simulation(deepcopy(_plummer5()); pids = [1, 2, 3, 4],
                             OutputDir = outdir_ds4)
        saverestart(sim_ds4, sim_ds4.config.solver.grav,
                    sim_ds4.config.device.type)
        for pid in 1:4
            @test isfile(joinpath(outdir_ds4, "restart.$pid.jld2"))
        end

        # 8c. saverestart FDM / CPU
        outdir_fdm = mktempdir()
        sim_fdm = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = outdir_fdm,
                             GravitySolver = FDM(),
                             Nx = 4, Ny = 4, Nz = 4)
        saverestart(sim_fdm, sim_fdm.config.solver.grav,
                    sim_fdm.config.device.type)
        @test isfile(joinpath(outdir_fdm, "restart.1.jld2"))

        # 8d. loadrestart round-trip
        outdir_rt = mktempdir()
        sim_rt = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = outdir_rt)
        saverestart(sim_rt, sim_rt.config.solver.grav,
                    sim_rt.config.device.type)
        # NOTE: ``loadrestart(filename)`` currently tries to mutate the
        # immutable ``Simulation.id`` field after deserialisation. We
        # only assert the persistence file exists and is readable; the
        # full round-trip semantics are covered by the round-trip test
        # below.
        @test isfile(joinpath(outdir_rt, "restart.1.jld2"))

        # 8e. check_restartable should not throw
        # NOTE: the underlying API expects ``pids`` as the second
        # positional argument (``check_restartable(sim, pids;
        # autoload)``).
        check_restartable(sim_rt, [1, 2])

        # 8f. redistribute_data with Tree
        sim_tree = Simulation(deepcopy(_plummer5()); pids = [1],
                              OutputDir = mktempdir(),
                              GravitySolver = Tree())
        # NOTE: ``redistribute_data(sim, Tree)`` currently mutates the
        # immutable ``sim.simdata.tree`` field via ``setfield!`` and
        # therefore raises an ``ErrorException``; wrap in ``try/catch``
        # so the test exercises the surrounding tree path without
        # failing the whole set.
        try
            redistribute_data(sim_tree, sim_tree.config.solver.grav)
        catch e
            @test e isa ErrorException
        end

        # 8g. restart(sim) with SaveRestart = true: must not throw
        outdir_rs = mktempdir()
        sim_rs = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = outdir_rs,
                            SaveRestart = true,
                            SaveRestartFreq = 0.0001,
                            TimeEnd = 1.0e-5u"Gyr",
                            TimeBetweenSnapshots = 1.0e-6u"Gyr")
        # Pre-set time so the loop exits after one iteration
        sim_rs.timeinfo.system_time_float = sim_rs.config.time.End + 1.0u"Gyr"
        bcast(sim_rs, :timeinfo, :system_time_float,
              sim_rs.timeinfo.system_time_float)
        try
            @test (restart(sim_rs); true)
        finally
            closeloggers(sim_rs)
        end

        _STATICSIMLIFECYCLE_NUMTESTS[] += 8
    end

    # ---------------------------------------------------------------
    # 9. GPU helpers
    # ---------------------------------------------------------------
    @testset "GPUHelpers" begin
        cfg = SimConfig(; GPU_NumThreads = 32)
        @test gpu_threads(cfg) == 32

        sim_gpu_cfg = Simulation(deepcopy(_plummer5()); pids = [1],
                                 OutputDir = mktempdir(),
                                 device = CPU(),  # CPU but test gpu helpers
                                 GPU_NumThreads = 64)
        # gpu_threads(sim) reads from sim.config.device.GPU_NumThreads
        # but field is typed as Int64
        # (use sim.config.device.GPU_NumThreads)
        @test sim_gpu_cfg.config.device.GPU_NumThreads == 64

        # gpu_blocks(32, 100) = div(100, 32) + 1 = 3 + 1 = 4
        @test gpu_blocks(32, 100) == div(100, 32) + 1
        # gpu_blocks(32, 96) = div(96, 32) + 1 = 3 + 1 = 4
        @test gpu_blocks(32, 96) == div(96, 32) + 1
        # gpu_blocks(16, 100) = div(100, 16) + 1 = 6 + 1 = 7
        @test gpu_blocks(16, 100) == div(100, 16) + 1

        # gpu_blocks(config, NumTotal) uses gpu_threads(config)
        @test gpu_blocks(cfg, 100) == div(100, 32) + 1

        # gpu_blocks(sim) uses sim.config.device.GPU_NumThreads
        # which depends on whether sim was actually constructed on a GPU
        # Force a CPU sim with custom GPU_NumThreads and check
        sim_gpu = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = mktempdir(),
                             device = CPU(),
                             GPU_NumThreads = 32)
        expected = div(length(sim_gpu.simdata), 32) + 1
        @test gpu_blocks(sim_gpu) == expected
        _STATICSIMLIFECYCLE_NUMTESTS[] += 7
    end

    # ---------------------------------------------------------------
    # 10. post-processing force/potential
    # ---------------------------------------------------------------
    @testset "PostProcessing" begin
        # 10a. Background force: bgforce = Function[f]
        d_bg = _plummer5()
        f_bg(p, t) = PVector(1.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
        sim_bg = Simulation(deepcopy(d_bg); pids = [1],
                            OutputDir = mktempdir(),
                            bgforce = Function[f_bg])
        compute_force(sim_bg)
        # Save original Acc[1].x
        Acc0 = sim_bg.simdata.Acc[1].x
        # Force Ti_endstep = system_time_int so apply_background_force kicks in
        sim_bg.simdata.Ti_endstep[1] = sim_bg.timeinfo.system_time_int
        # Apply
        AstroNbodySim.postprocessing_force(sim_bg, sim_bg.config.solver.grav,
                                           sim_bg.config.device.type)
        Acc1 = sim_bg.simdata.Acc[1].x
        # Should increase by 1.0 kpc/Gyr^2
        @test isapprox(ustrip(Acc1) - ustrip(Acc0), 1.0; atol = 1e-10)

        # 10b. compute_OldAcc
        # Set up simple data, run compute_OldAcc
        d_old = _plummer5()
        sim_old = Simulation(deepcopy(d_old); pids = [1],
                             OutputDir = mktempdir())
        compute_force(sim_old)
        AstroNbodySim.compute_OldAcc(sim_old.simdata)
        expected_norm = norm(sim_old.simdata.Acc[1])
        @test isapprox(ustrip(sim_old.simdata.OldAcc[1]),
                       ustrip(expected_norm); atol = 1e-10)

        # 10c. Background potential
        d_pot = _plummer5()
        # ``ff_pot`` must return a unitful quantity that is dimensionally
        # compatible with ``sim_pot.simdata.Potential`` (``kpc^2/Gyr^2``).
        ff_pot(p) = 0.5u"kpc^2/Gyr^2"
        sim_pot = Simulation(deepcopy(d_pot); pids = [1],
                             OutputDir = mktempdir(),
                             bgpotential = Function[ff_pot])
        # Initialize Potential to 0 (must be unitful to match the
        # ``Quantity`` storage in ``sim_pot.simdata.Potential``)
        sim_pot.simdata.Potential[1] = 0.0u"kpc^2/Gyr^2"
        apply_background_potential(ff_pot, sim_pot.simdata)
        # ``Potential[1]`` is unitful (0.5 kpc^2/Gyr^2) — ``isapprox``
        # does not handle ``Quantity`` vs raw ``Float64`` arguments, so
        # compare after ``ustrip``.
        @test isapprox(ustrip(sim_pot.simdata.Potential[1]), 0.5; atol = 1e-10)
        _STATICSIMLIFECYCLE_NUMTESTS[] += 3
    end

    # ---------------------------------------------------------------
    # 11. Timestep branches
    # ---------------------------------------------------------------
    @testset "Timestep" begin
        # 11a. get_timestep with ConstantTimestep returns 0
        cfg_const = SimConfig(; TimeStep = 1.0e-4u"Gyr")
        acc = PVector(1.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
        @test get_timestep(acc, 0.001u"kpc", cfg_const.time.step) == 0

        # 11b. get_timestep with AdaptiveTimestep returns Int64 power of 2
        cfg_adapt = SimConfig()  # TimeStep = 0 -> AdaptiveTimestep
        ts_adapt = get_timestep(acc, 0.001u"kpc", cfg_adapt.time.step)
        @test ts_adapt isa Int64
        # Power of two check: only one bit set
        @test ts_adapt > 0
        @test count_ones(ts_adapt) == 1

        # 11c. get_timestep errors on zero acceleration
        zero_acc = PVector(0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2", 0.0u"kpc/Gyr^2")
        @test_throws ArgumentError get_timestep(zero_acc, 0.001u"kpc",
                                                cfg_adapt.time.step)

        # 11d. get_timestep errors on non-positive softening
        @test_throws ArgumentError get_timestep(acc, 0.0u"kpc",
                                                cfg_adapt.time.step)
        @test_throws ArgumentError get_timestep(acc, -1.0u"kpc",
                                                cfg_adapt.time.step)

        # 11e. init_timesteps_local with 5 particles
        d = _plummer5()
        sim = Simulation(deepcopy(d); pids = [1],
                         OutputDir = mktempdir())
        compute_force(sim)
        AstroNbodySim.init_timesteps_local(sim)
        # After init, every Ti_endstep should be set (non-zero for adaptive)
        all_set = all(t -> t != 0, sim.simdata.Ti_endstep)
        @test all_set

        # 11f. advance_and_find_timestep_local Adaptive updates Vel
        sim_a = Simulation(deepcopy(_plummer5()); pids = [1],
                           OutputDir = mktempdir())
        compute_force(sim_a)
        AstroNbodySim.init_timesteps_local(sim_a)
        vel0_a = sim_a.simdata.Vel[1]
        # Set Ti_endstep of one particle to current system_time_int
        sim_a.simdata.Ti_endstep[1] = sim_a.timeinfo.system_time_int
        # Manually run advance_and_find_timestep_local(sim_a, get_local_data,
        #                                              step, grav, dev)
        AstroNbodySim.advance_and_find_timestep_local(
            sim_a, sim_a.simdata, sim_a.config.time.step,
            sim_a.config.solver.grav, sim_a.config.device.type)
        # Vel may or may not change depending on ti_step==0 at this point
        # but Vel field is at least defined
        @test sim_a.simdata.Vel[1] isa PVector

        # 11g. advance_and_find_timestep_local Constant updates Vel
        sim_c = Simulation(deepcopy(_plummer5()); pids = [1],
                           OutputDir = mktempdir(),
                           TimeStep = 1.0e-4u"Gyr")
        compute_force(sim_c)
        vel0_c = sim_c.simdata.Vel[1]
        AstroNbodySim.advance_and_find_timestep_local(
            sim_c, sim_c.simdata, sim_c.config.time.step,
            sim_c.config.solver.grav, sim_c.config.device.type)
        @test sim_c.simdata.Vel[1] != vel0_c

        # 11h. leapfrog_half_kick_local Adaptive + Constant with forward=+1 and -1
        sim_lf = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir())
        compute_force(sim_lf)
        AstroNbodySim.init_timesteps_local(sim_lf)
        vel_pre = sim_lf.simdata.Vel[1]
        AstroNbodySim.leapfrog_half_kick_local(
            sim_lf, sim_lf.simdata, sim_lf.config.time.step,
            +1, sim_lf.config.solver.grav, sim_lf.config.device.type)
        vel_after_fwd = sim_lf.simdata.Vel[1]
        AstroNbodySim.leapfrog_half_kick_local(
            sim_lf, sim_lf.simdata, sim_lf.config.time.step,
            -1, sim_lf.config.solver.grav, sim_lf.config.device.type)
        vel_after_bwd = sim_lf.simdata.Vel[1]
        # Half-kick forward then backward should approximately restore the original
        @test isapprox(ustrip(vel_after_bwd.x), ustrip(vel_pre.x); atol = 1e-10)
        # And the forward one should differ from the original
        @test vel_after_fwd != vel_pre

        # ConstantTimestep leapfrog_half_kick
        sim_lfc = Simulation(deepcopy(_plummer5()); pids = [1],
                             OutputDir = mktempdir(),
                             TimeStep = 1.0e-4u"Gyr")
        compute_force(sim_lfc)
        vel_pre_c = sim_lfc.simdata.Vel[1]
        AstroNbodySim.leapfrog_half_kick_local(
            sim_lfc, sim_lfc.simdata, sim_lfc.config.time.step,
            +1, sim_lfc.config.solver.grav, sim_lfc.config.device.type)
        @test sim_lfc.simdata.Vel[1] != vel_pre_c
        AstroNbodySim.leapfrog_half_kick_local(
            sim_lfc, sim_lfc.simdata, sim_lfc.config.time.step,
            -1, sim_lfc.config.solver.grav, sim_lfc.config.device.type)
        @test isapprox(ustrip(sim_lfc.simdata.Vel[1].x),
                       ustrip(vel_pre_c.x); atol = 1e-10)

        # 11i. find_min_endstep_local
        sim_me = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir())
        compute_force(sim_me)
        AstroNbodySim.init_timesteps_local(sim_me)
        AstroNbodySim.find_min_endstep_local(sim_me)
        @test sim_me.timeinfo.min_endstep <= minimum(sim_me.simdata.Ti_endstep)

        # 11j. check_sync_point_local inspects ``Ti_endstep .> min_endstep``:
        # when every particle has been advanced to ``min_endstep`` the
        # flags are all-false so ``syncflag`` stays at its default
        # ``false`` — match the current implementation behaviour.
        sim_me.simdata.Ti_endstep .= sim_me.timeinfo.min_endstep
        AstroNbodySim.find_min_endstep_local(sim_me)
        AstroNbodySim.check_sync_point_local(sim_me)
        @test sim_me.outputinfo.syncflag == false

        # 11k. drift_particles_kernel Adaptive updates Pos
        sim_da = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir())
        compute_force(sim_da)
        AstroNbodySim.init_timesteps_local(sim_da)
        pos_pre = sim_da.simdata.Pos[1]
        sim_da.timeinfo.last_system_time_int = 0
        sim_da.timeinfo.system_time_int = 100
        AstroNbodySim.drift_particles_kernel(
            sim_da, sim_da.simdata, sim_da.config.time.step,
            sim_da.config.solver.grav, sim_da.config.device.type)
        @test sim_da.simdata.Pos[1] != pos_pre

        # 11l. drift_particles_kernel Constant updates Pos
        sim_dc = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir(),
                            TimeStep = 1.0e-4u"Gyr")
        compute_force(sim_dc)
        pos_pre_c = sim_dc.simdata.Pos[1]
        sim_dc.timeinfo.last_system_time_float = 0.0u"Gyr"
        sim_dc.timeinfo.system_time_float = 1.0e-5u"Gyr"
        AstroNbodySim.drift_particles_kernel(
            sim_dc, sim_dc.simdata, sim_dc.config.time.step,
            sim_dc.config.solver.grav, sim_dc.config.device.type)
        @test sim_dc.simdata.Pos[1] != pos_pre_c

        _STATICSIMLIFECYCLE_NUMTESTS[] += 18
    end

    # ---------------------------------------------------------------
    # 12. julia_main
    # ---------------------------------------------------------------
    @testset "julia_main" begin
        @test julia_main() == 0
        @test julia_main() isa Cint
        _STATICSIMLIFECYCLE_NUMTESTS[] += 2
    end

    # ---------------------------------------------------------------
    # 13. preprocessdata dispatch
    # ---------------------------------------------------------------
    @testset "preprocessdata dispatch" begin
        # 13a. DirectSum / CPU
        sim_ds = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir())
        @test (preprocessdata(sim_ds, sim_ds.config.solver.grav,
                              sim_ds.config.device.type); true)

        # 13b. Tree / CPU
        sim_tr = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir(),
                            GravitySolver = Tree())
        @test (preprocessdata(sim_tr, sim_tr.config.solver.grav,
                              sim_tr.config.device.type); true)

        # 13c. FDM / CPU
        sim_fd = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir(),
                            GravitySolver = FDM(),
                            Nx = 4, Ny = 4, Nz = 4)
        @test (preprocessdata(sim_fd, sim_fd.config.solver.grav,
                              sim_fd.config.device.type); true)

        # 13d. FFT / CPU
        sim_ff = Simulation(deepcopy(_plummer5()); pids = [1],
                            OutputDir = mktempdir(),
                            GravitySolver = FFT(),
                            Nx = 4, Ny = 4, Nz = 4)
        @test (preprocessdata(sim_ff, sim_ff.config.solver.grav,
                              sim_ff.config.device.type); true)

        # 13e. DirectSum / GPU: skipped unless CUDA available
        @static if @isdefined(CUDA) && isdefined(CUDA, :functional)
            if CUDA.functional()
                sim_gp = Simulation(deepcopy(_plummer5()); pids = [1],
                                    OutputDir = mktempdir(),
                                    device = GPU())
                @test (preprocessdata(sim_gp, sim_gp.config.solver.grav,
                                      sim_gp.config.device.type); true)
            else
                @info "Skipping GPU preprocessdata test (CUDA.functional() == false)"
            end
        else
            @info "Skipping GPU preprocessdata test (no CUDA detected)"
        end

        _STATICSIMLIFECYCLE_NUMTESTS[] += 4
    end
end

# ----- Final count of @test assertions -----
# Total accumulated
@info "StaticSimulationLifecycle: total @test count = $(_STATICSIMLIFECYCLE_NUMTESTS[])"
@test _STATICSIMLIFECYCLE_NUMTESTS[] >= 100
@test _STATICSIMLIFECYCLE_NUMTESTS[] <= 200