@setup_workload begin
    @compile_workload begin
        binary_rot_vel(mass, D) = sqrt(G * mass / D / 2.0)
        binary_orbit_period(mass, D) = 2 * pi * sqrt(0.5 * D^3 / G / mass)

        function compute_dt(a::Number; SofteningLength = 0.0001u"kpc", ErrTolTimestep = 0.025)
            return sqrt(2 * ErrTolTimestep * SofteningLength / abs(a))
        end

        p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
        p1y(sim::Simulation) = find_particle(sim, 1).Pos.y
        analysers = Dict(
            "x" => p1x,
            "y" => p1y,
        )
        
        G = Constant(uAstro).G

        mass = 1.0e10u"Msun"
        D = 2.0u"kpc"
        vel = binary_rot_vel(mass, D)
        T = binary_orbit_period(mass, D)

        data = StructArray(Star(uAstro, id = i) for i in 1:2);
        data.Pos[1] = PVector(+0.5D, 0.0u"kpc", 0.0u"kpc");
        data.Pos[2] = PVector(-0.5D, 0.0u"kpc", 0.0u"kpc");
        data.Vel[1] = PVector(0.0u"kpc/Gyr", +vel, 0.0u"kpc/Gyr");
        data.Vel[2] = PVector(0.0u"kpc/Gyr", -vel, 0.0u"kpc/Gyr");
        data.Mass .= mass;

        TimeStep = compute_dt(G*mass/D^2)
        TimeEnd = T * 0.01
        TimeBetweenSnapshots = TimeEnd * 0.5

        ds = Simulation(
            deepcopy(data);
            analysers,
            TimeEnd,
            TimeStep,
            TimeBetweenSnapshots,
            OutputDir = joinpath(tempdir(), "AstroNbodySim.jl"),
        );
        run(ds)
        compute_potential(ds)

        suggest_softlen!(ds)

        # tree = Simulation(
        #     deepcopy(data);
        #     analysers,
        #     TimeEnd,
        #     TimeStep,
        #     TimeBetweenSnapshots,
        #     GravitySolver = Tree(),
        #     OutputDir = joinpath(tempdir(), "AstroNbodySim.jl"),
        # );
        # run(tree)
    end
end
