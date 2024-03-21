mkpathIfNotExist("output/PM")

@testset "Particle Mesh" begin
    function plot_phi_slice(sim::Simulation, units)
        L = size(sim.simdata.phi)[3]
        U = isnothing(units) ? "Unitless" : "Unitful"
        S = sim.config.solver.grav isa FFT ? 2 : 1
        for i in S:L
            p = plot_slice(sim.simdata, :phi, i, units; title = "Potential slice")
            Makie.save("output/PM/SingleParticle-$(traitstring(sim.config.solver.grav))-$(U)-PotentialSlice-$(i).png", p)
        end

        return true
    end

    @testset "Single-particle test" begin
        d = StructArray([Star(nothing)])
        d.Mass[1] = 1.0e10

        ds = Simulation(deepcopy(d);
            units = nothing,
            constants = Constant(nothing, uAstro),
        )

        N = 10

        @testset "FDM" begin
            sim = Simulation(deepcopy(d);
                units = nothing,
                constants = Constant(nothing, uAstro),
                GravitySolver = FDM(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                xMin = -2.0,
                xMax = +2.0,
                yMin = -2.0,
                yMax = +2.0,
                zMin = -2.0,
                zMax = +2.0,
                NG = 1,
            )

            compute_force(sim)
            @test plot_phi_slice(sim, nothing)
            
            # Mass conservation
            @test sum(sim.simdata.rho) * prod(sim.simdata.config.Δ) ≈ 1.0e10

            center_N = div(N,2) + 2
            center_rho = sim.simdata.rho[center_N,center_N,center_N]
            center_mass = center_rho * prod(sim.simdata.config.Δ)
            @test center_mass ≈ 1.0e10

            # Potential
            G = sim.config.constants.G
            h = 0.0001
            Pos = sim.simdata.pos[2,2,2]
            
            ## compare to DirectSum
            Pot = compute_potential(ds, Pos, h)[1]
            @test abs(sim.simdata.phi[2,2,2] - Pot) / abs(Pot) < 0.001

            ## compare to direct sum in mesh
            center_pos = sim.simdata.pos[center_N,center_N,center_N]
            Pot = -G * center_mass / sqrt((Pos-center_pos)*(Pos-center_pos) + h^2)
            @test abs(sim.simdata.phi[2,2,2] - Pot) / abs(Pot) < 0.001
        end

        @testset "FFT" begin
            sim = Simulation(deepcopy(d);
                units = nothing,
                constants = Constant(nothing, uAstro),
                GravitySolver = FFT(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                xMin = -2.0,
                xMax = +2.0,
                yMin = -2.0,
                yMax = +2.0,
                zMin = -2.0,
                zMax = +2.0,
                NG = 1,
                BoundaryCondition = Dirichlet(),
            )

            compute_force(sim)
            @test plot_phi_slice(sim, nothing)
        end
    end

    @testset "Single-particle test (unitful)" begin
        d = StructArray([Star(uAstro)])
        d.Mass[1] = 1.0e10u"Msun"

        ds = Simulation(deepcopy(d);)

        N = 10

        @testset "FDM" begin
            sim = Simulation(deepcopy(d);
                GravitySolver = FDM(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                xMin = -2.0u"kpc",
                xMax = +2.0u"kpc",
                yMin = -2.0u"kpc",
                yMax = +2.0u"kpc",
                zMin = -2.0u"kpc",
                zMax = +2.0u"kpc",
                NG = 1,
            )
            compute_force(sim)
            @test plot_phi_slice(sim, uAstro)
            
            @test sum(sim.simdata.rho) * prod(sim.simdata.config.Δ) ≈ 1.0e10u"Msun"

            center_N = div(N,2) + 2
            center_rho = sim.simdata.rho[center_N,center_N,center_N]
            center_mass = center_rho * prod(sim.simdata.config.Δ)
            @test center_mass ≈ 1.0e10u"Msun"

            # Potential
            G = sim.config.constants.G
            h = 0.0001u"kpc"
            Pos = sim.simdata.pos[2,2,2]
            
            ## compare to DirectSum
            Pot = compute_potential(ds, Pos, h)[1]
            @test abs(sim.simdata.phi[2,2,2] - Pot) / abs(Pot) < 0.001

            ## compare to direct sum in mesh
            center_pos = sim.simdata.pos[center_N,center_N,center_N]
            Pot = -G * center_mass / sqrt((Pos-center_pos)*(Pos-center_pos) + h^2)
            @test abs(sim.simdata.phi[2,2,2] - Pot) / abs(Pot) < 0.001
        end

        @testset "FFT" begin
            sim = Simulation(deepcopy(d);
                GravitySolver = FFT(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                xMin = -2.0u"kpc",
                xMax = +2.0u"kpc",
                yMin = -2.0u"kpc",
                yMax = +2.0u"kpc",
                zMin = -2.0u"kpc",
                zMax = +2.0u"kpc",
                NG = 1,
                BoundaryCondition = Dirichlet(),
            )
            compute_force(sim)
            @test plot_phi_slice(sim, uAstro)
        end
    end

    @testset "Force curve" begin
        R = collect(0.0:0.01:2.0)
        pos = [PVector(x, 0.0, 0.0, u"kpc") for x in R]
        
        N = 10
        @testset "Single-particle" begin
            d = StructArray([Star(uAstro)])
            d.Mass[1] = 1.0e10u"Msun"

            ds = Simulation(deepcopy(d);)

            pm = Simulation(deepcopy(d);
                GravitySolver = FDM(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                xMin = -2.0u"kpc",
                xMax = +2.0u"kpc",
                yMin = -2.0u"kpc",
                yMax = +2.0u"kpc",
                zMin = -2.0u"kpc",
                zMax = +2.0u"kpc",
                NG = 1,
            )
            compute_force(pm)

            acc_ds = norm.(compute_force(ds, pos, 0.01u"kpc"))
            acc_pm = norm.(compute_force(pm, pos, 0.01u"kpc"))
            
            f = Figure(; resolution = (800,450))
            axis = GLMakie.Axis(f[1,1],
                title = "Radial acceleration of single particle",
                xlabel = "log10(R [kpc])",
                ylabel = "log10(acc [kpc/Gyr^2])",
            )
            p1 = GLMakie.lines!(axis, log10.(R), log10.(ustrip.(u"kpc/Gyr^2", acc_ds)))
            p2 = GLMakie.lines!(axis, log10.(R), log10.(ustrip.(u"kpc/Gyr^2", acc_pm)))
            legend = GLMakie.Legend(
                f[1,1], [p1,p2], ["DirectSum", "Particle-Mesh"],
                tellheight = false,
                tellwidth = false,
                halign = :right,
                valign = :top,
                margin = (10, 10, 10, 10),
            )
            Makie.save("output/PM/StaticForceSingleParticle.png", f)
        end

        @testset "Plummer star cluster" begin
            d = generate(PlummerStarCluster(
                NumSamples = 1000,
                VirialRadius = 0.5u"kpc"
            ))
            
            ds = Simulation(deepcopy(d);)
            compute_force(ds)

            pm = Simulation(deepcopy(d);
                GravitySolver = FDM(),
                device = CPU(),
                Nx = N,
                Ny = N,
                Nz = N,
                EnlargeMesh = 2.0,
                NG = 1,
            )
            compute_force(pm)


            f = Figure(; resolution = (800,450))
            axis = GLMakie.Axis(f[1,1],
                title = "Radial acceleration of Plummer star cluster",
                xlabel = "log10(R [kpc])",
                ylabel = "log10(acc [kpc/Gyr^2])",
            )
            p1 = plot_force!(axis, get_local_data(ds), uAstro; savelog=false, markersize = 3.0)
            p2 = plot_force!(axis, get_local_data(pm), uAstro; savelog=false, markersize = 3.0)
            legend = GLMakie.Legend(
                f[1,1], [p1,p2], ["DirectSum", "Particle-Mesh"],
                tellheight = false,
                tellwidth = false,
                halign = :right,
                valign = :top,
                margin = (10, 10, 10, 10),
            )
            Makie.save("output/PM/StaticForcePlummer.png", f)
        end
    end
end

sol(x::Real) = sin(2*pi*x) + sin(32*pi*x) / 256
init_rho(x::Real) = -4*pi*pi*sin(2*pi*x) - 4*pi*pi*sin(32*pi*x)

sol(p::PVector2D) =  sin(2*pi*p.x) * sin(2*pi*p.y) + sin(32*pi*p.x) * sin(32*pi*p.y) / 256
init_rho(p::PVector2D) = -8 * pi * pi * sin(2*pi*p.x) * sin(2*pi*p.y) - 8 * pi * pi * sin(32*pi*p.x) * sin(32*pi*p.y)

sol(p::PVector) =  sin(2*pi*p.x) * sin(2*pi*p.y) * sin(2*pi*p.z) + sin(32*pi*p.x) * sin(32*pi*p.y) * sin(2*pi*p.z) / 256
init_rho(p::PVector) = -12 * pi * pi * sin(2*pi*p.x) * sin(2*pi*p.y) * sin(2*pi*p.z) - 12 * pi * pi * sin(32*pi*p.x) * sin(32*pi*p.y) * sin(32*pi*p.z)

L2norm(r) = sqrt(sum((r.^2)/prod(size(r))))

@testset "fft poisson solver" begin
    @testset "1D" begin
        function test_fft1D(Nx; boundary = Periodic())
            m = MeshCartesianStatic(;
                xMin = 0.0,
                xMax = 1.0,
                Nx = Nx - 1,
                NG = 0,
                dim = 1,
                boundary,
            )
            m.rho .= init_rho.(m.pos)
            fft_poisson!(m, m.rho, m.config.boundary)
            s = sol.(m.pos)
            r = m.phi .- s

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Error of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, r)
            Makie.save("output/PM/fft1D-$(traitstring(boundary))-$(Nx)-error.png", f)

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Solution of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, m.phi)
            Makie.save("output/PM/fft1D-$(traitstring(boundary))-$(Nx)-solution.png", f)

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Exact value of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, s)
            Makie.save("output/PM/fft1D-$(traitstring(boundary))-$(Nx)-exact.png", f)

            return L2norm(r)
        end

        #! FFT mesh should be at least larger than 32x32

        @test test_fft1D(8; boundary = Periodic()) < 0.5
        @test test_fft1D(16; boundary = Periodic()) < 0.89
        @test test_fft1D(32; boundary = Periodic()) < 0.06
        @test test_fft1D(64; boundary = Periodic()) < 0.031
        @test test_fft1D(128; boundary = Periodic()) < 0.016
        @test test_fft1D(256; boundary = Periodic()) < 0.008
        @test test_fft1D(512; boundary = Periodic()) < 0.004
        @test test_fft1D(1024; boundary = Periodic()) < 0.002
        @test test_fft1D(2048; boundary = Periodic()) < 0.001
        
        
        @test test_fft1D(8; boundary = Dirichlet()) < 0.7
        @test test_fft1D(16; boundary = Dirichlet()) < 1.0
        @test test_fft1D(32; boundary = Dirichlet()) < 0.12
        @test test_fft1D(64; boundary = Dirichlet()) < 0.06
        @test test_fft1D(128; boundary = Dirichlet()) < 0.031
        @test test_fft1D(256; boundary = Dirichlet()) < 0.016
        @test test_fft1D(512; boundary = Dirichlet()) < 0.008
        @test test_fft1D(1024; boundary = Dirichlet()) < 0.004
        @test test_fft1D(2048; boundary = Dirichlet()) < 0.002
    end

    @testset "2D" begin
        function test_fft2D(Nx, Ny = Nx; boundary = Periodic())
            m = MeshCartesianStatic(;
                xMin = 0.0,
                yMin = 0.0,
                xMax = 1.0,
                yMax = 1.0,
                Nx = Nx - 1,
                Ny = Ny - 1,
                NG = 0,
                dim = 2,
                boundary,
            )
            m.rho .= init_rho.(m.pos)
            fft_poisson!(m, m.rho, m.config.boundary)
            s = sol.(m.pos)
            r = m.phi .- s

            f = plot_mesh_heatmap(r, title = "Error of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fft2D-$(traitstring(boundary))-$(Nx)x$(Ny)-error.png", f)

            f = plot_mesh_heatmap(m.phi, title = "Solution of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fft2D-$(traitstring(boundary))-$(Nx)x$(Ny)-solution.png", f)

            f = plot_mesh_heatmap(s, title = "Exact value of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fft2D-$(traitstring(boundary))-$(Nx)x$(Ny)-exact.png", f)
            return L2norm(r)
        end
        
        #! FFT mesh should be at least larger than 32x32

        @test test_fft2D(8; boundary = Periodic()) < 0.3
        @test test_fft2D(16; boundary = Periodic()) < 0.61
        @test test_fft2D(32; boundary = Periodic()) < 0.05
        @test test_fft2D(64; boundary = Periodic()) < 0.03
        @test test_fft2D(128; boundary = Periodic()) < 0.012
        @test test_fft2D(256; boundary = Periodic()) < 0.006
        @test test_fft2D(512; boundary = Periodic()) < 0.003
        @test test_fft2D(1024; boundary = Periodic()) < 0.002
        @test test_fft2D(2048; boundary = Periodic()) < 0.0007
        
        
        @test test_fft2D(8; boundary = Dirichlet()) < 0.3
        @test test_fft2D(16; boundary = Dirichlet()) < 0.61
        @test test_fft2D(32; boundary = Dirichlet()) < 0.06
        @test test_fft2D(64; boundary = Dirichlet()) < 0.03
        @test test_fft2D(128; boundary = Dirichlet()) < 0.02
        @test test_fft2D(256; boundary = Dirichlet()) < 0.01
        @test test_fft2D(512; boundary = Dirichlet()) < 0.004
        @test test_fft2D(1024; boundary = Dirichlet()) < 0.002
        @test test_fft2D(2048; boundary = Dirichlet()) < 0.001
    end

    @testset "3D" begin
        Nx = Ny = Nz = 256
        boundary = Periodic()
        
        function test_fft3D(Nx, Ny = Nx, Nz = Nx; boundary = Periodic())
            m = MeshCartesianStatic(;
                xMin = 0.0,
                yMin = 0.0,
                zMin = 0.0,
                xMax = 1.0,
                yMax = 1.0,
                zMax = 1.0,
                Nx = Nx - 1,
                Ny = Ny - 1,
                Nz = Nz - 1,
                NG = 0,
                dim = 3,
                boundary,
            )
            m.rho .= init_rho.(m.pos)
            fft_poisson!(m, m.rho, m.config.boundary)
            s = sol.(m.pos)
            r = m.phi .- s

            f = GLMakie.volume(r)
            Makie.save("output/PM/fft3D-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-error.png", f)

            f = GLMakie.volume(m.phi)
            Makie.save("output/PM/fft3D-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-solution.png", f)

            f = GLMakie.volume(s)
            Makie.save("output/PM/fft3D-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-exact.png", f)

            return L2norm(r)
        end

        #! FFT mesh should be at least larger than 32x32

        @test test_fft3D(8; boundary = Periodic()) < 0.15
        @test test_fft3D(16; boundary = Periodic()) < 0.41
        @test test_fft3D(32; boundary = Periodic()) < 0.04
        @test test_fft3D(64; boundary = Periodic()) < 0.02
        @test test_fft3D(128; boundary = Periodic()) < 0.01
        @test test_fft3D(256; boundary = Periodic()) < 0.005
        # It's becoming slow and taking too much memory
        
        @test test_fft3D(8; boundary = Dirichlet()) < 0.14
        @test test_fft3D(16; boundary = Dirichlet()) < 0.40
        @test test_fft3D(32; boundary = Dirichlet()) < 0.04
        @test test_fft3D(64; boundary = Dirichlet()) < 0.021
        @test test_fft3D(128; boundary = Dirichlet()) < 0.011
        @test test_fft3D(256; boundary = Dirichlet()) < 0.01
        # It's becoming slow and taking too much memory
    end
end

@testset "fdm poisson solver" begin
    @testset "1D" begin
        function test_fdm1D(Nx; boundary = Periodic(), device = CPU(), sparse = true)
            m = MeshCartesianStatic(;
                xMin = 0.0,
                xMax = 1.0,
                Nx = Nx - 1,
                NG = 0,
                dim = 1,
                boundary,
                device,
            )
            m.rho .= init_rho.(m.pos)
            fdm_poisson(m, Val(1), 1/(4*pi), m.config.mode, device, boundary, sparse)
            s = sol.(Array(m.pos))
            r = Array(m.phi) .- s

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Error of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, r)
            Makie.save("output/PM/fdm1D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)-error.png", f)

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Solution of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, Array(m.phi))
            Makie.save("output/PM/fdm1D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)-solution.png", f)

            f = Figure(; resolution = (1200,900))
            a = GLMakie.Axis(f[1,1]; title = "Exact value of potential $(Nx) $(traitstring(boundary))")
            GLMakie.lines!(a, s)
            Makie.save("output/PM/fdm1D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)-exact.png", f)

            return L2norm(r)
        end

        @test test_fdm1D(8; boundary = Periodic(), device = CPU(), sparse = true) < 0.92
        @test test_fdm1D(16; boundary = Periodic(), device = CPU(), sparse = true) < 0.91
        @test test_fdm1D(32; boundary = Periodic(), device = CPU(), sparse = true) < 0.1
        @test test_fdm1D(64; boundary = Periodic(), device = CPU(), sparse = true) < 0.04
        @test test_fdm1D(128; boundary = Periodic(), device = CPU(), sparse = true) < 0.02
        @test test_fdm1D(256; boundary = Periodic(), device = CPU(), sparse = true) < 0.013
        @test test_fdm1D(512; boundary = Periodic(), device = CPU(), sparse = true) < 0.005
        @test test_fdm1D(1024; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(2048; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(4096; boundary = Periodic(), device = CPU(), sparse = true) < 0.0006
        @test test_fdm1D(8192; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(16384; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(32768; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(65536; boundary = Periodic(), device = CPU(), sparse = true) < 0.0006
        
        # It's becoming slow
        #@test test_fdm1D(131072; boundary = Periodic(), device = CPU(), sparse = true) < 0.002
        #@test test_fdm1D(262144; boundary = Periodic(), device = CPU(), sparse = true) < 0.002

        @test test_fdm1D(8; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.7
        @test test_fdm1D(16; boundary = Dirichlet(), device = CPU(), sparse = true) < 1.07
        @test test_fdm1D(32; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.12
        @test test_fdm1D(64; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.06
        @test test_fdm1D(128; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.04
        @test test_fdm1D(256; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.016
        @test test_fdm1D(512; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.008
        @test test_fdm1D(1024; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.004
        @test test_fdm1D(2048; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.002
        @test test_fdm1D(4096; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.001
        @test test_fdm1D(8192; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0005
        @test test_fdm1D(16384; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0003
        @test test_fdm1D(32768; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0002
        @test test_fdm1D(65536; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        @test test_fdm1D(131072; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        @test test_fdm1D(262144; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        @test test_fdm1D(524288; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        @test test_fdm1D(1048576; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        @test test_fdm1D(2097152; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0001
        
        # It's becoming slow and taking too much memory
        #@test test_fdm1D(4194304; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0002
        #@test test_fdm1D(8388608; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.0006

        # We test GPU and dense array to make sure that it is functional
        @test test_fdm1D(8; boundary = Periodic(), device = GPU(), sparse = false) < 6.4
        @test test_fdm1D(16; boundary = Periodic(), device = GPU(), sparse = false) < 9.2
        @test test_fdm1D(32; boundary = Periodic(), device = GPU(), sparse = false) < 181
        @test test_fdm1D(64; boundary = Periodic(), device = GPU(), sparse = false) < 0.43
        @test test_fdm1D(128; boundary = Periodic(), device = GPU(), sparse = false) < 0.81
        @test test_fdm1D(256; boundary = Periodic(), device = GPU(), sparse = false) < 0.93
        @test test_fdm1D(512; boundary = Periodic(), device = GPU(), sparse = false) < 1.4
        @test test_fdm1D(1024; boundary = Periodic(), device = GPU(), sparse = false) < 913
        @test test_fdm1D(2048; boundary = Periodic(), device = GPU(), sparse = false) < 0.7
        @test test_fdm1D(4096; boundary = Periodic(), device = GPU(), sparse = false) < 0.43
        @test test_fdm1D(8192; boundary = Periodic(), device = GPU(), sparse = false) < 0.5

        # It's becoming slow and taking too much memory
        #@test test_fdm1D(16384; boundary = Periodic(), device = GPU(), sparse = false) < 0.3
        #@test test_fdm1D(32768; boundary = Periodic(), device = GPU(), sparse = false) < 0.3

        @test test_fdm1D(8; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.7
        @test test_fdm1D(16; boundary = Dirichlet(), device = GPU(), sparse = false) < 1.07
        @test test_fdm1D(32; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.12
        @test test_fdm1D(64; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.06
        @test test_fdm1D(128; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.04
        @test test_fdm1D(256; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.016
        @test test_fdm1D(512; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.008
        @test test_fdm1D(1024; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.004
        @test test_fdm1D(2048; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.002
        @test test_fdm1D(4096; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.005
        @test test_fdm1D(8192; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.12

        # Vacuum boundary condition is tested in Particle-Mesh simulations
    end

    @testset "2D" begin
        function test_fdm2D(Nx, Ny = Nx; boundary = Periodic(), device = CPU(), sparse = true)
            m = MeshCartesianStatic(;
                xMin = 0.0,
                yMin = 0.0,
                xMax = 1.0,
                yMax = 1.0,
                Nx = Nx - 1,
                Ny = Ny - 1,
                NG = 0,
                dim = 2,
                boundary,
                device,
            )
            m.rho .= init_rho.(m.pos)
            fdm_poisson(m, Val(2), 1/(4*pi), m.config.mode, device, boundary, sparse)
            s = CUDA.@allowscalar sol.(Array(m.pos))
            r = CUDA.@allowscalar Array(m.phi) .- s

            f = plot_mesh_heatmap(r, title = "Error of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fdm2D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)-error.png", f)

            f = CUDA.@allowscalar plot_mesh_heatmap(Array(m.phi), title = "Solution of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fdm2D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)-solution.png", f)

            f = plot_mesh_heatmap(s, title = "Exact value of potential $(Nx)x$(Ny) $(traitstring(boundary))")
            Makie.save("output/PM/fdm2D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)-exact.png", f)

            return L2norm(r)
        end

        @test test_fdm2D(8; boundary = Periodic(), device = CPU(), sparse = true) < 0.24
        @test test_fdm2D(16; boundary = Periodic(), device = CPU(), sparse = true) < 0.61
        @test test_fdm2D(32; boundary = Periodic(), device = CPU(), sparse = true) < 0.05
        @test test_fdm2D(64; boundary = Periodic(), device = CPU(), sparse = true) < 0.03
        @test test_fdm2D(128; boundary = Periodic(), device = CPU(), sparse = true) < 0.02
        @test test_fdm2D(256; boundary = Periodic(), device = CPU(), sparse = true) < 0.006
        @test test_fdm2D(512; boundary = Periodic(), device = CPU(), sparse = true) < 0.003
        # It's becoming slow and taking too much memory
        #@test test_fdm2D(1024; boundary = Periodic(), device = CPU(), sparse = true) < 0.003
        
        @test test_fdm2D(8; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.33
        @test test_fdm2D(16; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.71
        @test test_fdm2D(32; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.09
        @test test_fdm2D(64; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.05
        @test test_fdm2D(128; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.03
        @test test_fdm2D(256; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.02
        @test test_fdm2D(512; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.006
        # It's becoming slow and taking too much memory

        # We test GPU and dense array to make sure that it is functional
        @test test_fdm2D(8; boundary = Periodic(), device = GPU(), sparse = false) < 1.43
        @test test_fdm2D(16; boundary = Periodic(), device = GPU(), sparse = false) < 0.97
        @test test_fdm2D(32; boundary = Periodic(), device = GPU(), sparse = false) < 0.22
        @test test_fdm2D(64; boundary = Periodic(), device = GPU(), sparse = false) < 0.07
        @test test_fdm2D(128; boundary = Periodic(), device = GPU(), sparse = false) < 0.03

        @test test_fdm2D(8; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.33
        @test test_fdm2D(16; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.71
        @test test_fdm2D(32; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.09
        @test test_fdm2D(64; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.05
        @test test_fdm2D(128; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.03

        # Vacuum boundary condition is tested in Particle-Mesh simulations
    end

    @testset "3D" begin
        function test_fdm3D(Nx, Ny = Nx, Nz = Nx; boundary = Periodic(), device = CPU(), sparse = true)
            m = MeshCartesianStatic(;
                xMin = 0.0,
                yMin = 0.0,
                zMin = 0.0,
                xMax = 1.0,
                yMax = 1.0,
                zMax = 1.0,
                Nx = Nx - 1,
                Ny = Ny - 1,
                Nz = Nz - 1,
                NG = 0,
                dim = 3,
                boundary,
                device,
            )
            m.rho .= init_rho.(m.pos)
            fdm_poisson(m, Val(3), 1/(4*pi), m.config.mode, device, boundary, sparse)
            s = CUDA.@allowscalar sol.(Array(m.pos))
            r = CUDA.@allowscalar Array(m.phi) .- s

            f = GLMakie.volume(r)
            Makie.save("output/PM/fdm3D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-error.png", f)

            f = CUDA.@allowscalar GLMakie.volume(Array(m.phi))
            Makie.save("output/PM/fdm3D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-solution.png", f)

            f = GLMakie.volume(s)
            Makie.save("output/PM/fdm3D-$(traitstring(device))-$(traitstring(boundary))-$(Nx)x$(Ny)x$(Nz)-exact.png", f)

            return L2norm(r)
        end

        @test test_fdm3D(8; boundary = Periodic(), device = CPU(), sparse = true) < 0.15
        @test test_fdm3D(16; boundary = Periodic(), device = CPU(), sparse = true) < 0.41
        @test test_fdm3D(24; boundary = Periodic(), device = CPU(), sparse = true) < 0.11
        # It's becoming slow
        #@test test_fdm3D(32; boundary = Periodic(), device = CPU(), sparse = true) < 1.54 #! why larger?
        
        @test test_fdm3D(8; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.19
        @test test_fdm3D(16; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.48
        @test test_fdm3D(24; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.08
        # It's becoming slow
        #@test test_fdm3D(32; boundary = Dirichlet(), device = CPU(), sparse = true) < 0.061

        # We test GPU and dense array to make sure that it is functional
        @test test_fdm3D(8; boundary = Periodic(), device = GPU(), sparse = false) < 0.16
        @test test_fdm3D(16; boundary = Periodic(), device = GPU(), sparse = false) < 0.41
        @test test_fdm3D(18; boundary = Periodic(), device = GPU(), sparse = false) < 0.33
        @test test_fdm3D(20; boundary = Periodic(), device = GPU(), sparse = false) < 0.07
        @test test_fdm3D(22; boundary = Periodic(), device = GPU(), sparse = false) < 0.051

        @test test_fdm3D(8; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.19
        @test test_fdm3D(16; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.48
        @test test_fdm3D(18; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.33
        @test test_fdm3D(20; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.1
        @test test_fdm3D(22; boundary = Dirichlet(), device = GPU(), sparse = false) < 0.1

        # Vacuum boundary condition is tested in Particle-Mesh simulations
    end
end