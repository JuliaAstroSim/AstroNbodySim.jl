#! This is a heavy-load test. Run at local computers!
mkpathIfNotExist("output/Timestep/")

# when r > h, SofteningLength is not working

function compute_acc(x::Number; SofteningLength = 0.0001u"kpc", G = Constant(uAstro).G, M = 1.0e10u"Msun")
    return -G * M / (4 * x^2#= + SofteningLength^2=#)
end

function compute_dt(a::Number; SofteningLength = 0.0001u"kpc", ErrTolTimestep = 0.025)
    return sqrt(2 * ErrTolTimestep * SofteningLength / abs(a))
end

function compute_acc(x::PVector; SofteningLength = 0.0001u"kpc", G = Constant(uAstro).G, M = 1.0e10u"Msun")
    return -G * M / (4 * x * x#= + SofteningLength^2=#) * normalize(ustrip(x))
end

function compute_dt(a::PVector; SofteningLength = 0.0001u"kpc", ErrTolTimestep = 0.025)
    return sqrt(2 * ErrTolTimestep * SofteningLength / norm(a))
end

function integrate_const_euler(N, t1, x1, v1;
        compute_acc::Function = compute_acc,
        compute_dt::Function = compute_dt,
        filename = "output/Timestep/ConstEuler.csv",
    )
    a1 = compute_acc(x1)
    dt = compute_dt(a1)
    df = DataFrame(
        t = Vector{typeof(t1)}(undef, N),
        x = Vector{typeof(x1)}(undef, N),
        v = Vector{typeof(v1)}(undef, N),
        a = Vector{typeof(a1)}(undef, N),
        dt = Vector{typeof(dt)}(undef, N),
    )
    df.t[1] = t1
    df.x[1] = x1
    df.v[1] = v1
    df.a[1] = a1
    df.dt[1] = dt

    for i in 1:N-1
        df.x[i+1] = df.x[i] + df.v[i] * dt
        df.v[i+1] = df.v[i] + df.a[i] * dt
        df.a[i+1] = compute_acc(df.x[i+1])
        df.t[i+1] = df.t[i] + dt
        df.dt[i+1] = dt
    end

    CSV.write(filename, df)
    return df
end

function integrate_const_leapfrog(N, t1, x1, v1;
        compute_acc::Function = compute_acc,
        compute_dt::Function = compute_dt,
        filename = "output/Timestep/ConstLeapfrog.csv",
    )
    a1 = compute_acc(x1)
    dt = compute_dt(a1)
    df = DataFrame(
        t = Vector{typeof(t1)}(undef, N),
        x = Vector{typeof(x1)}(undef, N),
        v = Vector{typeof(v1)}(undef, N),
        a = Vector{typeof(a1)}(undef, N),
        dt = Vector{typeof(dt)}(undef, N),
    )
    df.t[1] = t1
    df.x[1] = x1
    df.v[1] = v1
    df.a[1] = a1
    df.dt[1] = dt

    # forward half kick: v_{1 + 1/2} = v_1 0.5 a * dt. We rename v_1 to v_{1 + 1/2}
    df.v[1] += 0.5 * df.a[1] * dt

    # same as euler
    for i in 1:N-1
        df.x[i+1] = df.x[i] + df.v[i] * dt
        df.a[i+1] = compute_acc(df.x[i+1])
        df.v[i+1] = df.v[i] + df.a[i+1] * dt
        df.t[i+1] = df.t[i] + dt
        df.dt[i+1] = dt
    end

    # backward half kick: v_{n} here is actually v_{n + 1/2}. We kick it a half step backward
    df.v[end] -= 0.5 * df.a[end] * dt

    CSV.write(filename, df)
    return df
end

function integrate_adapt_euler(N, t1, x1, v1;
        compute_acc::Function = compute_acc,
        compute_dt::Function = compute_dt,
        filename = "output/Timestep/AdaptEuler.csv",
    )
    a1 = compute_acc(x1)
    dt = compute_dt(a1)
    df = DataFrame(
        t = Vector{typeof(t1)}(undef, N),
        x = Vector{typeof(x1)}(undef, N),
        v = Vector{typeof(v1)}(undef, N),
        a = Vector{typeof(a1)}(undef, N),
        dt = Vector{typeof(dt)}(undef, N),
    )
    df.t[1] = t1
    df.x[1] = x1
    df.v[1] = v1
    df.a[1] = a1
    df.dt[1] = dt

    for i in 1:N-1
        df.x[i+1] = df.x[i] + df.v[i] * df.dt[i]
        df.v[i+1] = df.v[i] + df.a[i] * df.dt[i]
        df.a[i+1] = compute_acc(df.x[i+1])
        df.t[i+1] = df.t[i] + df.dt[i]      # Adapt
        df.dt[i+1] = compute_dt(df.a[i+1])  # Adapt
    end

    CSV.write(filename, df)
    return df
end

function integrate_adapt_leapfrog(N, t1, x1, v1;
        compute_acc::Function = compute_acc,
        compute_dt::Function = compute_dt,
        filename = "output/Timestep/AdaptEuler.csv",
    )
    a1 = compute_acc(x1)
    dt = compute_dt(a1)
    df = DataFrame(
        t = Vector{typeof(t1)}(undef, N),
        x = Vector{typeof(x1)}(undef, N),
        v = Vector{typeof(v1)}(undef, N),
        a = Vector{typeof(a1)}(undef, N),
        dt = Vector{typeof(dt)}(undef, N),
    )
    df.t[1] = t1
    df.x[1] = x1
    df.v[1] = v1
    df.a[1] = a1
    df.dt[1] = dt

    # forward half kick: v_{1 + 1/2} = v_1 0.5 a * dt. We rename v_1 to v_{1 + 1/2}
    df.v[1] += 0.5 * df.a[1] * df.dt[1]

    for i in 1:N-1
        df.x[i+1] = df.x[i] + df.v[i] * df.dt[i]
        df.a[i+1] = compute_acc(df.x[i+1])
        df.dt[i+1] = compute_dt(df.a[i+1])
        df.v[i+1] = df.v[i] + df.a[i+1] * 0.5 * (df.dt[i] + df.dt[i+1])
        df.t[i+1] = df.t[i] + df.dt[i]
    end

    # backward half kick: v_{n} here is actually v_{n + 1/2}. We kick it a half step backward
    df.v[end] -= 0.5 * df.a[end] * df.dt[end]

    CSV.write(filename, df)
    return df
end

@everywhere p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
@everywhere p1y(sim::Simulation) = find_particle(sim, 1).Pos.y

function plot_orbit(sim, title::String)
    df = DataFrame(CSV.File(joinpath(sim.config.output.dir, "analysis.csv")))
    scene, layout = layoutscene(resolution = (960, 1080))
    axis = layout[1,1] = GLMakie.Axis(scene; title, xlabel = "x [kpc]", ylabel = "y [kpc]")
    axis.autolimitaspect = 1
    Makie.lines!(axis, df.x, df.y)
    Makie.save(joinpath("output/Timestep", title * ".png"), scene);
    return df
end

@testset "Time integration: binary free fall" begin
    G = Constant(uAstro).G
    M = 1.0e10u"Msun"
    R = 1.0u"kpc"

    ErrTolTimestep = 0.025
    SofteningLength = 0.0001u"kpc"

    d = StructArray(Star(uAstro, id = i) for i in 1:2);
    d.Mass .= M
    d.Pos[1] = PVector(+R, 0.0u"kpc", 0.0u"kpc")
    d.Pos[2] = PVector(-R, 0.0u"kpc", 0.0u"kpc")

    t1 = 0.0u"Gyr"
    x1 = +R
    v1 = -100.0u"kpc/Gyr"

    N = 260
    
    # Exact solutions
       ConstEuler = integrate_const_euler(   N, t1, x1, v1, filename = "output/Timestep/BinaryFreeFall-ConstEuler.csv")
    ConstLeapfrog = integrate_const_leapfrog(N, t1, x1, v1, filename = "output/Timestep/BinaryFreeFall-ConstLeapfrog.csv")
       AdaptEuler = integrate_adapt_euler(   N, t1, x1, v1, filename = "output/Timestep/BinaryFreeFall-AdaptEuler.csv")
    AdaptLeapfrog = integrate_adapt_leapfrog(N, t1, x1, v1, filename = "output/Timestep/BinaryFreeFall-AdaptLeapfrog.csv")


    f = CairoMakie.Figure(; resolution = (1600, 900))
    ax = CairoMakie.Axis(f[1,1],
        title = "Time integration test: binary free fall",
        xlabel = "t [Gyr]",
        ylabel = "x [kpc]",
    )
    CE = CairoMakie.lines!(ax, ustrip.(u"Gyr",    ConstEuler.t), ustrip.(u"kpc", ConstEuler.x),    color = colors[1])
    CL = CairoMakie.lines!(ax, ustrip.(u"Gyr", ConstLeapfrog.t), ustrip.(u"kpc", ConstLeapfrog.x), color = colors[2])
    AE = CairoMakie.lines!(ax, ustrip.(u"Gyr",    AdaptEuler.t), ustrip.(u"kpc", AdaptEuler.x),    color = colors[3])
    AL = CairoMakie.lines!(ax, ustrip.(u"Gyr", AdaptLeapfrog.t), ustrip.(u"kpc", AdaptLeapfrog.x), color = colors[4])
    legend = CairoMakie.Legend(f[1,2],
        [CE, CL, AE, AL],
        ["Euler const", "Leapfrog const", "Euler adapt", "Leapfrog adapt"],
    )
    Makie.save("output/Timestep/TimeIntegration-BinaryFreeFall.png", f)
end

@testset "Time integration: circular binary" begin
    G = Constant(uAstro).G
    M = 1.0e10u"Msun"
    D = 2.0u"kpc"

    ErrTolTimestep = 0.025
    SofteningLength = 0.0001u"kpc"

    binary_orbit_period(M, D) = 2 * pi * sqrt(0.5 * D^3 / G / M)

    binary_rot_vel(M, D) = sqrt(G * M / D / 2.0)

    vel = binary_rot_vel(M, D)
    T = binary_orbit_period(M, D)

    # Define the two particles
    d = StructArray(Star(uAstro, id = i) for i in 1:2);
    d.Pos[1] = PVector(+0.5D, 0.0u"kpc", 0.0u"kpc");
    d.Pos[2] = PVector(-0.5D, 0.0u"kpc", 0.0u"kpc");
    d.Vel[1] = PVector(0.0u"kpc/Gyr", +vel, 0.0u"kpc/Gyr");
    d.Vel[2] = PVector(0.0u"kpc/Gyr", -vel, 0.0u"kpc/Gyr");
    d.Mass .= M;

    t1 = 0.0u"Gyr"
    x1 = d.Pos[1]
    v1 = d.Vel[1]
    a1 = compute_acc(x1)

    N = 2900

    # Exact solutions
       ConstEuler = integrate_const_euler(   N, t1, x1, v1, filename = "output/Timestep/BinaryCircular-ConstEuler.csv")
    ConstLeapfrog = integrate_const_leapfrog(N, t1, x1, v1, filename = "output/Timestep/BinaryCircular-ConstLeapfrog.csv")
       AdaptEuler = integrate_adapt_euler(   N, t1, x1, v1, filename = "output/Timestep/BinaryCircular-AdaptEuler.csv")
    AdaptLeapfrog = integrate_adapt_leapfrog(N, t1, x1, v1, filename = "output/Timestep/BinaryCircular-AdaptLeapfrog.csv")

    f = CairoMakie.Figure(; resolution = (1000, 1000))
    ax = CairoMakie.Axis(f[1,1],
        title = "Time integration test: binary circular orbit (constant timestep)",
        xlabel = "x [kpc]",
        ylabel = "y [kpc]",
        aspect = AxisAspect(1.0),
    )
    CE = CairoMakie.lines!(ax, ustrip.(u"kpc", StructArray(ConstEuler.x).x),    ustrip.(u"kpc", StructArray(   ConstEuler.x).y), color = colors[1],)
    CL = CairoMakie.lines!(ax, ustrip.(u"kpc", StructArray(ConstLeapfrog.x).x), ustrip.(u"kpc", StructArray(ConstLeapfrog.x).y), color = colors[2],)
    CairoMakie.xlims!(ax, (-1.1,1.1))
    CairoMakie.ylims!(ax, (-1.1,1.1))
    legend = CairoMakie.Legend(f[1,2],
        [CE, CL],
        ["Euler const", "Leapfrog const"],
    )
    Makie.save("output/Timestep/TimeIntegration-BinaryCircular-Const.png", f)

    f = CairoMakie.Figure(; resolution = (1000, 1000))
    ax = CairoMakie.Axis(f[1,1],
        title = "Time integration test: binary circular orbit (adaptive timestep)",
        xlabel = "x [kpc]",
        ylabel = "y [kpc]",
        aspect = AxisAspect(1.0),
    )
    AE = CairoMakie.lines!(ax, ustrip.(u"kpc", StructArray(AdaptEuler.x).x),    ustrip.(u"kpc", StructArray(   AdaptEuler.x).y), color = colors[3],)
    AL = CairoMakie.lines!(ax, ustrip.(u"kpc", StructArray(AdaptLeapfrog.x).x), ustrip.(u"kpc", StructArray(AdaptLeapfrog.x).y), color = colors[4],)
    CairoMakie.xlims!(ax, (-1.1,1.1))
    CairoMakie.ylims!(ax, (-1.1,1.1))
    legend = CairoMakie.Legend(f[1,2],
        [AE, AL],
        ["Euler adapt", "Leapfrog adapt"],
    )
    Makie.save("output/Timestep/TimeIntegration-BinaryCircular-Adapt.png", f)

    analysers = Dict(
        "x" => p1x,
        "y" => p1y,
    )

    @testset "Constant Euler" begin
        TimeEnd = ConstEuler.t[end]
        TimeBetweenSnapshots = TimeEnd
        TimeStep = compute_dt(a1)
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-DirectSum-Serial")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-DirectSum-MAW")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-DirectSum-MTW")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-Tree-Serial")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-Tree-MAW")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-ConstEuler-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstEuler-Tree-MTW")
            @test df.x[end] ≈ ConstEuler.x[end].x
            @test df.y[end] ≈ ConstEuler.x[end].y
        end
    end

    @testset "Constant Leapfrog" begin
        TimeEnd = ConstLeapfrog.t[end]
        TimeBetweenSnapshots = TimeEnd
        TimeStep = compute_dt(a1)
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-DirectSum-Serial")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-DirectSum-MAW")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-DirectSum-MTW")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-Tree-Serial")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-Tree-MAW")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-ConstLeapfrog-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-ConstLeapfrog-Tree-MTW")
            @test df.x[end] ≈ ConstLeapfrog.x[end].x
            @test df.y[end] ≈ ConstLeapfrog.x[end].y
        end
    end

    @testset "Adaptive Euler" begin
        TimeEnd = AdaptEuler.t[end]
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-DirectSum-Serial")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-DirectSum-MAW")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-DirectSum-MTW")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-Tree-Serial")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-Tree-MAW")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptEuler-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptEuler-Tree-MTW")
            @test df.x[end] ≈ AdaptEuler.x[end].x
            @test df.y[end] ≈ AdaptEuler.x[end].y
        end
    end


    @testset "Adaptive Leapfrog" begin
        TimeEnd = AdaptLeapfrog.t[end]
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-DirectSum-Serial")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-DirectSum-MAW")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-DirectSum-MTW")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-Tree-Serial")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-Tree-MAW")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircular-AdaptLeapfrog-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircular-AdaptLeapfrog-Tree-MTW")
            @test df.x[end] ≈ AdaptLeapfrog.x[end].x
            @test df.y[end] ≈ AdaptLeapfrog.x[end].y
        end
    end
end


@testset "Time integration: circular binary background force" begin
    @everywhere G = Constant(uAstro).G
    @everywhere mass = 0.5e10u"Msun"
    R = 1.0u"kpc"

    ErrTolTimestep = 0.025
    SofteningLength = 0.0001u"kpc"

    binary_orbit_period(mass, R) = 2 * pi * sqrt(R^3 / G / mass)
    binary_rot_vel(mass, R) = sqrt(G * mass / R)

    vel = binary_rot_vel(mass, R)
    T = binary_orbit_period(mass, R)

    println("Estimated Period: ", T)
    println("Estimated Rotation Vel: ", vel)

    # Define the background force field
    @everywhere attractor(p::AbstractParticle) = -G * mass / (p.Pos * p.Pos) * normalize(p.Pos) / 1.0u"kpc"


    # Define two massless particles
    d = StructArray(Star(uAstro, id = i) for i in 1:2);
    d.Pos[1] = PVector(+R, 0.0u"kpc", 0.0u"kpc");
    d.Pos[2] = PVector(-R, 0.0u"kpc", 0.0u"kpc");
    d.Vel[1] = PVector(0.0u"kpc/Gyr", +vel, 0.0u"kpc/Gyr");
    d.Vel[2] = PVector(0.0u"kpc/Gyr", -vel, 0.0u"kpc/Gyr");

    bgforce = Function[attractor]
    analysers = Dict(
        "x" => p1x,
        "y" => p1y,
    )

    @testset "Constant Euler" begin
        TimeEnd = T
        TimeStep = compute_dt(attractor(d[1]))
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-DirectSum-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-DirectSum-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-DirectSum-MTW")
            @test !isnothing(df)
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-Tree-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-Tree-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstEuler-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstEuler-Tree-MTW")
            @test !isnothing(df)
        end
    end

    @testset "Constant Leapfrog" begin
        TimeEnd = T
        TimeStep = compute_dt(attractor(d[1]))
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-DirectSum-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-DirectSum-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-DirectSum-MTW")
            @test !isnothing(df)
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-Tree-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-Tree-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep,
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-ConstLeapfrog-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-ConstLeapfrog-Tree-MTW")
            @test !isnothing(df)
        end
    end

    @testset "Adaptive Euler" begin
        TimeEnd = T
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-DirectSum-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-DirectSum-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-DirectSum-MTW")
            @test !isnothing(df)
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-Tree-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-Tree-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Euler(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptEuler-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptEuler-Tree-MTW")
            @test !isnothing(df)
        end
    end


    @testset "Adaptive Leapfrog" begin
        TimeEnd = T
        TimeBetweenSnapshots = TimeEnd
        @testset "DirectSum" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-DirectSum-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-DirectSum-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-DirectSum-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-DirectSum-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-DirectSum-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-DirectSum-MTW")
            @test !isnothing(df)
        end

        @testset "Tree" begin
            sim = Simulation(
                deepcopy(d),
                pids = [1];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-Tree-Serial",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-Tree-Serial")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [1,2];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-Tree-MAW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-Tree-MAW")
            @test !isnothing(df)

            sim = Simulation(
                deepcopy(d),
                pids = [2,3];
                analysers,
                bgforce,
                TimeEnd,
                TimeStep = 0.0u"Gyr",
                TimeBetweenSnapshots,
                GravitySolver = Tree(),
                TimeIntegrationAlgorithm = Leapfrog(),
                OutputDir = "output/Timestep/BinaryCircularBG-AdaptLeapfrog-Tree-MTW",
            )
            set_softlen!(sim, 0.0001u"kpc")
            run(sim)
            df = plot_orbit(sim, "BinaryCircularBG-AdaptLeapfrog-Tree-MTW")
            @test !isnothing(df)
        end
    end
end