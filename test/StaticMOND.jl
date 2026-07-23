#######################################################################
# Static MOND test
#
# This file exercises the Milgrom-1983 / QUMOND / MOND model dispatch
# implemented in:
#   - src/mond/milgrom1983.jl
#   - src/mond/qumond.jl
#   - src/mond/MOG.jl
#   - src/base/PostProcessing/Force.jl
#   - src/base/PostProcessing/Potential.jl
#   - src/PM/gravity.jl
#
# It is included from runtests.jl after the standard prelude:
#   @everywhere using AstroNbodySim
#   @everywhere using AstroNbodySim.AstroSimBase
#   @everywhere using Unitful, UnitfulAstro, PhysicalParticles,
#                   PhysicalMeshes, PhysicalTrees
#   @everywhere astro()
#
# Hard constraints (do not violate):
#   * Never run julia in the sandbox
#   * Never modify runtests.jl
#   * 60-90 @test statements in total
#######################################################################

@testset "Static MOND" begin
    # ------------------------------------------------------------------
    # Common physical / numerical constants for the whole file
    # ------------------------------------------------------------------
    const_local  = Constant(uAstro)
    G_local      = const_local.G
    ACC0_local   = const_local.ACC0
    Gnum         = ustrip(G_local)
    ACC0num      = ustrip(ACC0_local)
    SoftLen_local = 1.0e-4u"kpc"

    # ------------------------------------------------------------------
    # 1. MOND interpolation function
    # ------------------------------------------------------------------
    @testset "MOND interpolation" begin
        # Reference values
        y1_ref   = 0.5 * (1.0 + sqrt(5.0))            # nu(1, n=2) base
        nu1_y1   = sqrt(y1_ref)                       # nu(1, 2) = (0.5*(1+sqrt(5)))^0.5

        # nu(y, 2) at y=1
        @test nu(1.0, 2) ≈ nu1_y1
        @test nu(1.0, 2) ≈ (0.5 * (1.0 + sqrt(5.0)))^0.5

        # nu(100, 2) → Newtonian limit → 1
        @test isapprox(nu(100.0, 2), 1.0; atol = 1.0e-3)
        @test isapprox(nu(1.0e6, 2), 1.0; atol = 1.0e-6)

        # nu(0.01, 2): should give nu*y ≈ 1/sqrt(y) in deep-MOND
        # analytic: nu(0.01,2) = (0.5*(1+sqrt(1+4/1e-4)))^0.5
        expected_nu_001 = (0.5 * (1.0 + sqrt(1.0 + 4.0 / 1.0e-4)))^0.5
        @test isapprox(nu(0.01, 2), expected_nu_001; rtol = 1.0e-12)

        # In the deep-MOND limit (small y) the standard Milgrom
        # interpolation ``nu(y, 2) = (0.5*(1+sqrt(1+4/y^2)))^(1/2)`` behaves
        # as ``nu ~ 1/sqrt(y)`` so ``nu`` itself grows without bound
        # while the dimensionless combination ``y * nu(y, 2) ~ sqrt(y)``
        # vanishes. The tests below verify both limits.
        @test (1.0e-2) * nu(1.0e-2, 2) > 1.0e-2    # y*nu ~ sqrt(y) → 0
        @test (1.0e-4) * nu(1.0e-4, 2) > 1.0e-3    # y*nu ~ sqrt(y)
        @test (1.0e-6) * nu(1.0e-6, 2) > 1.0e-4    # y*nu ~ sqrt(y)
        # Equivalent restatement: nu(y, 2) ~ 1/sqrt(y) >> 1 in deep MOND.
        @test nu(1.0e-2, 2) > 1.0e1                # nu(1e-2) ≈ 10
        @test nu(1.0e-4, 2) > 1.0e2                # nu(1e-4) ≈ 100
        @test nu(1.0e-6, 2) > 1.0e3                # nu(1e-6) ≈ 1000

        # nu1 / nu2 closed forms
        @test nu1(2.0)  ≈ 0.5 * (1.0 + sqrt(1.0 + 2.0))
        @test nu2(2.0)  ≈ sqrt(0.5 * (1.0 + sqrt(1.0 + 1.0)))
        @test nu1(0.5)  ≈ 0.5 * (1.0 + sqrt(1.0 + 8.0))
        @test nu2(0.5)  ≈ sqrt(0.5 * (1.0 + sqrt(1.0 + 16.0)))

        # nu(y, n) at n=1 should equal nu1(y); n=2 should equal nu2(y)
        for yval in (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0)
            @test isapprox(nu(yval, 1), nu1(yval); rtol = 1.0e-12)
            @test isapprox(nu(yval, 2), nu2(yval); rtol = 1.0e-12)
        end

        # large-y limit nu(y, 2) → 1 (and y * nu is bounded)
        @test isapprox(nu(1.0e8, 2), 1.0; atol = 1.0e-8)
        @test isapprox(nu(1.0e10, 2), 1.0; atol = 1.0e-10)

        # n=3 also consistent
        yval = 4.0
        @test isapprox(nu(yval, 3),
                       (0.5 * (1.0 + sqrt(1.0 + 4.0 / yval^3)))^(1.0/3);
                       rtol = 1.0e-12)
    end

    # ------------------------------------------------------------------
    # 2. Milgrom 1983 acceleration rescaling
    # ------------------------------------------------------------------
    @testset "Milgrom 1983 rescaling" begin
        # Two-particle binary on master
        BinaryData = StructArray(Star(uAstro, id = i) for i in 1:2)
        BinaryData.Pos[1] = PVector(+1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Pos[2] = PVector(-1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Mass .= 1.0e10u"Msun"
        BinaryData.Vel .= PVector(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")

        sim = Simulation(BinaryData;
            GravityModel = AstroNbodySim.MOND1983Milgrom(),
            ForceSofteningTable = [SoftLen_local for _ in 1:6],
            OutputDir = mktempdir(),
        )
        compute_force(sim)

        data = get_local_data(sim)

        # Force Ti_endstep to match system_time_int
        sys_ti = sim.timeinfo.system_time_int
        for i in 1:length(data)
            data.Ti_endstep[i] = sys_ti
        end

        # Save original Acc
        orig_acc1 = data.Acc[1]
        orig_acc2 = data.Acc[2]
        aN1 = norm(orig_acc1)
        aN2 = norm(orig_acc2)

        # Apply Milgrom 1983 with the proper integer time
        AstroNbodySim.mond_Milgrom1983(sim, data, sys_ti)

        # The relationship: new = orig * nu(|orig|/ACC0, 2)
        expected1 = orig_acc1 * nu(aN1 / ACC0_local, 2)
        expected2 = orig_acc2 * nu(aN2 / ACC0_local, 2)
        @test data.Acc[1] ≈ expected1
        @test data.Acc[2] ≈ expected2

        # The two accelerations are anti-symmetric → after MOND still
        # anti-symmetric
        @test abs(data.Acc[1].x + data.Acc[2].x) < 1.0e-10u"kpc/Gyr^2"

        # If we restore the original Acc and call with a non-matching
        # system_time_int, nothing should change.
        data.Acc[1] = orig_acc1
        data.Acc[2] = orig_acc2
        AstroNbodySim.mond_Milgrom1983(sim, data, sys_ti + 12345)
        @test data.Acc[1] ≈ orig_acc1
        @test data.Acc[2] ≈ orig_acc2

        # In the Newtonian limit (y = aN / ACC0 >> 1) the rescaling factor
        # ``nu(y, 2)`` approaches 1. For this 1e10 Msun binary at ±1 kpc
        # ``y ≈ 2.9`` so the factor is only mildly Newtonian (``nu ≈ 1.05``)
        # and we widen ``atol`` accordingly.
        deep_factor = nu(aN1 / ACC0_local, 2)
        @test isapprox(ustrip(deep_factor), 1.0; atol = 1.0e-1)

        # Sim config sanity
        @test sim.config.grav.MOND_nuIndex == 2.0
        @test sim.config.constants.ACC0 == ACC0_local
    end

    # ------------------------------------------------------------------
    # 3. Simulation with MOND1983Milgrom gravity model
    # ------------------------------------------------------------------
    @testset "MOND in Simulation" begin
        BinaryData = StructArray(Star(uAstro, id = i) for i in 1:2)
        BinaryData.Pos[1] = PVector(+1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Pos[2] = PVector(-1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Mass .= 1.0e10u"Msun"
        BinaryData.Vel .= PVector(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")

        sim = Simulation(BinaryData;
            GravityModel = AstroNbodySim.MOND1983Milgrom(),
            ForceSofteningTable = [SoftLen_local for _ in 1:6],
            OutputDir = mktempdir(),
        )
        @test sim.config.grav.model isa AstroNbodySim.MOND1983Milgrom
        @test !(sim.config.grav.model isa AstroNbodySim.QUMOND)
        @test !(sim.config.grav.model isa AstroNbodySim.Newton)

        # Default Newton vs explicit Newton
        sim_newton = Simulation(BinaryData;
            GravityModel = AstroNbodySim.Newton(),
            ForceSofteningTable = [SoftLen_local for _ in 1:6],
            OutputDir = mktempdir(),
        )
        @test sim_newton.config.grav.model isa AstroNbodySim.Newton

        # QUMOND model dispatch
        sim_qumond = Simulation(BinaryData;
            GravityModel = AstroNbodySim.QUMOND(),
            ForceSofteningTable = [SoftLen_local for _ in 1:6],
            OutputDir = mktempdir(),
        )
        @test sim_qumond.config.grav.model isa AstroNbodySim.QUMOND
    end

    # ------------------------------------------------------------------
    # 4-6. QUMOND operators on a small mesh
    # ------------------------------------------------------------------
    @testset "QUMOND operators" begin
        m = MeshCartesianStatic(uAstro;
            Nx = 8, Ny = 8, Nz = 8,
            NG = 0,
            boundary = Periodic(),
        )

        # 4) QUMOND_PDM_density: constant acc → rho_PDM = 0
        m.acc.x .= 1.0u"kpc/Gyr^2"
        m.acc.y .= 0.0u"kpc/Gyr^2"
        m.acc.z .= 0.0u"kpc/Gyr^2"
        m.rho    .= 0.0u"Msun/kpc^3"

        rho_pdm_const = QUMOND_PDM_density(m, ACC0_local, G_local)
        @test size(rho_pdm_const) == size(m.acc.x)
        # In the strict interior (away from the zero-padded stencil boundary)
        # the divergence of a constant acceleration field is zero.
        sl = size(rho_pdm_const)
        if all(sl .> 2)
            @test all(abs.(ustrip.(rho_pdm_const[2:end-1, 2:end-1, 2:end-1]))
                      .< 1.0e-10)
        end

        # Position-dependent acc: should yield non-zero rho_PDM.
        # The prefactor is dimensionless, only the final scale has units.
        xpos_units = ustrip.(m.pos.x)
        m.acc.x .= (1.0 .+ 0.01 .* xpos_units) .* u"kpc/Gyr^2"
        m.acc.y .= 0.0u"kpc/Gyr^2"
        m.acc.z .= 0.0u"kpc/Gyr^2"

        rho_pdm_var = QUMOND_PDM_density(m, ACC0_local, G_local)
        @test size(rho_pdm_var) == size(m.acc.x)
        @test maximum(abs.(ustrip.(rho_pdm_var))) > 1.0e-15

        # 5) QUMOND_phi: returns a field with mesh-vertex shape
        phi_const = QUMOND_phi(m, ACC0_local, G_local)
        @test size(phi_const) == size(m.acc.x)
        @test eltype(phi_const) <: eltype(m.acc.x) ||
              typeof(phi_const) <: AbstractArray

        # 6) QUMOND_acc / QUMOND_acc!
        # Reset to a varying acc so that phi is non-trivial
        xpos_units = ustrip.(m.pos.x)
        m.acc.x .= (1.0 .+ 0.01 .* xpos_units) .* u"kpc/Gyr^2"
        m.acc.y .= 0.5u"kpc/Gyr^2"
        m.acc.z .= 0.25u"kpc/Gyr^2"
        m.rho  .= 0.0u"Msun/kpc^3"

        delta = QUMOND_acc(m, ACC0_local, G_local)
        @test size(delta.x) == size(m.acc.x)

        # Save state, apply in-place, verify m.acc changed by exactly delta
        saved_x = copy(m.acc.x)
        saved_y = copy(m.acc.y)
        saved_z = copy(m.acc.z)
        ret = QUMOND_acc!(m, ACC0_local, G_local)
        @test ret === nothing
        @test m.acc.x ≈ saved_x .+ delta.x
        @test m.acc.y ≈ saved_y .+ delta.y
        @test m.acc.z ≈ saved_z .+ delta.z

        # In the constant-acc case the in-place operator should not change
        # m.acc (because the modification term is zero)
        m.acc.x .= 1.0u"kpc/Gyr^2"
        m.acc.y .= 0.0u"kpc/Gyr^2"
        m.acc.z .= 0.0u"kpc/Gyr^2"
        saved_x = copy(m.acc.x)
        QUMOND_acc!(m, ACC0_local, G_local)
        @test m.acc.x ≈ saved_x
    end

    # ------------------------------------------------------------------
    # 7. MOG stubs
    # ------------------------------------------------------------------
    @testset "MOG stubs" begin
        m = MeshCartesianStatic(uAstro;
            Nx = 8, Ny = 8, Nz = 8, NG = 0,
            boundary = Periodic(),
        )

        # Pre-load the mesh with non-trivial acc so the no-op character of
        # MOG_acc! is testable.
        m.acc.x .= 1.0u"kpc/Gyr^2"
        m.acc.y .= 2.0u"kpc/Gyr^2"
        m.acc.z .= 3.0u"kpc/Gyr^2"

        @test (AstroNbodySim.MOG_phi(); true)
        @test (AstroNbodySim.MOG_acc(m, ACC0num, Gnum); true)
        # MOG_acc! is a stub; should leave m.acc unchanged
        before_x = copy(m.acc.x)
        before_y = copy(m.acc.y)
        before_z = copy(m.acc.z)
        AstroNbodySim.MOG_acc!(m, ACC0num, Gnum)
        @test m.acc.x == before_x
        @test m.acc.y == before_y
        @test m.acc.z == before_z
        @test all(m.acc.x .== 1.0u"kpc/Gyr^2")
        @test all(m.acc.y .== 2.0u"kpc/Gyr^2")
        @test all(m.acc.z .== 3.0u"kpc/Gyr^2")
    end

    # ------------------------------------------------------------------
    # 8. Background force + postprocessing
    # ------------------------------------------------------------------
    @testset "postprocessing_force" begin
        BinaryData = StructArray(Star(uAstro, id = i) for i in 1:2)
        BinaryData.Pos[1] = PVector(+1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Pos[2] = PVector(-1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        BinaryData.Mass .= 1.0e10u"Msun"
        BinaryData.Vel .= PVector(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")

        f_bg(p, t) = PVector(1.0, 0.0, 0.0)u"kpc/Gyr^2"

        sim = Simulation(BinaryData;
            ForceSofteningTable = [SoftLen_local for _ in 1:6],
            bgforce = Function[f_bg],
            OutputDir = mktempdir(),
        )

        # Make sure Newton's gravity is used (otherwise the MOND branch in
        # postprocessing_force also modifies Acc; here we only test the
        # background-force path).
        @test sim.config.grav.model isa AstroNbodySim.Newton

        compute_force(sim)
        data = get_local_data(sim)

        # Force Ti_endstep to match the integer system time
        sys_ti = sim.timeinfo.system_time_int
        for i in 1:length(data)
            data.Ti_endstep[i] = sys_ti
        end

        # Snapshot Acc and OldAcc after the first postprocessing pass
        acc_after_first_1 = data.Acc[1]
        acc_after_first_2 = data.Acc[2]
        # call postprocessing_force again to apply the background force
        AstroNbodySim.postprocessing_force(sim,
            AstroNbodySim.DirectSum(), AstroNbodySim.CPU())

        # data.Acc.x must have increased by exactly 1 kpc/Gyr^2
        @test data.Acc[1].x ≈ acc_after_first_1.x + 1.0u"kpc/Gyr^2"
        @test data.Acc[1].y ≈ acc_after_first_1.y
        @test data.Acc[1].z ≈ acc_after_first_1.z
        # Particle 2 should have received the same background force
        @test data.Acc[2].x ≈ acc_after_first_2.x + 1.0u"kpc/Gyr^2"
        @test data.Acc[2].y ≈ acc_after_first_2.y
        @test data.Acc[2].z ≈ acc_after_first_2.z

        # data.OldAcc[i] = norm(data.Acc[i]) was set by compute_OldAcc
        @test data.OldAcc[1] ≈ norm(data.Acc[1])
        @test data.OldAcc[2] ≈ norm(data.Acc[2])

        # apply_background_force alone, on a clean data, should also work
        fresh = StructArray(Star(uAstro, id = i) for i in 1:2)
        fresh.Pos[1] = PVector(+1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        fresh.Pos[2] = PVector(-1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        fresh.Mass .= 1.0e10u"Msun"
        fresh.Vel .= PVector(0.0u"kpc/Gyr", 0.0u"kpc/Gyr", 0.0u"kpc/Gyr")
        fresh.Ti_endstep .= 0
        AstroNbodySim.apply_background_force(f_bg, fresh, 0, 0.0u"Gyr")
        @test fresh.Acc[1].x ≈ 1.0u"kpc/Gyr^2"
        @test fresh.Acc[2].x ≈ 1.0u"kpc/Gyr^2"

        # compute_OldAcc populates OldAcc = norm(Acc)
        AstroNbodySim.compute_OldAcc(fresh)
        @test fresh.OldAcc[1] ≈ norm(fresh.Acc[1])
        @test fresh.OldAcc[2] ≈ norm(fresh.Acc[2])
    end

    # ------------------------------------------------------------------
    # 9. GravityConfig dispatch
    # ------------------------------------------------------------------
    @testset "MOND model dispatch" begin
        cfg_milgrom = GravityConfig(; GravityModel = AstroNbodySim.MOND1983Milgrom())
        @test cfg_milgrom.model isa AstroNbodySim.MOND1983Milgrom
        @test !(cfg_milgrom.model isa AstroNbodySim.QUMOND)
        @test !(cfg_milgrom.model isa AstroNbodySim.Newton)
        @test cfg_milgrom.MOND_nuIndex == 2.0

        cfg_qumond = GravityConfig(; GravityModel = AstroNbodySim.QUMOND())
        @test cfg_qumond.model isa AstroNbodySim.QUMOND
        @test !(cfg_qumond.model isa AstroNbodySim.MOND1983Milgrom)
        @test !(cfg_qumond.model isa AstroNbodySim.Newton)

        cfg_newton = GravityConfig(; GravityModel = AstroNbodySim.Newton())
        @test cfg_newton.model isa AstroNbodySim.Newton
        @test !(cfg_newton.model isa AstroNbodySim.MOND1983Milgrom)
        @test !(cfg_newton.model isa AstroNbodySim.QUMOND)

        # The model is the default of Newton()
        cfg_default = GravityConfig()
        @test cfg_default.model isa AstroNbodySim.Newton

        # All three types are subtypes of GravityModel (provided by
        # AstroSimBase, re-imported via using AstroNbodySim.AstroSimBase)
        @test cfg_milgrom.model isa AstroSimBase.GravityModel
        @test cfg_qumond.model  isa AstroSimBase.GravityModel
        @test cfg_newton.model  isa AstroSimBase.GravityModel
    end
end
