# In this example, we lookback the trajectories of MW's satellites
# The original src/data/MW_satellites_origin.csv are extracted from battaglia2022gaia_Astron. Astrophys. - Gaia early DR3 systemic motions of local group dwarf galaxies and orbital properties with a massive Large Magellanic Cloud
#    and converted to Galactocentric coordinates using astropy.
# The parameters of LMC are from vasiliev2021tango_MNRAS - Tango for three sagittarius, LMC, and the milky way

using Unitful, UnitfulAstro
using StructArrays
using DataFrames, CSV
using GLMakie
using ProgressMeter

using AstroIC
using AstroNbodySim
using WaveDM

using PhysicalParticles


df_MW_satellites = AstroIC.load_data_MW_satellites()  # load from src/data/MW_satellites.csv


##### General parameters

Np = 1000
# Np = 5000
# Np = 10000
# Np = 50000
GravitySolver = Tree()
# GravitySolver = DirectSum()
pids = [1]

TimeEnd = 6.0u"Gyr"
TimeStep = 0.0u"Gyr"  # adaptive if zero.
# TimeStep = 1e-5u"Gyr"
TimeBetweenSnapshots = 0.1u"Gyr"

file_LMC = joinpath(@__DIR__, "output/LMC_traj_lookback_Vasiliev2021", "analysis.csv")
flag_load_LMC = false    # If false, re-compute the trajectory of LMC
# flag_load_LMC = true   # If true and the file exists, skip trajectory lookback of LMC
Np_LMC = 200
# Np_LMC = 1000
# Np_LMC = 10000

##### Initialize MW model
@info "Initializing MW's baryons"
#TODO: WaveDM.jl
baryon_particles = test_MW_MOND(;
    V = (x,y,z,ψ)->0.0,
    Nx = 32,
    Np,
    baryon_mode = :particles_static,
    export_particles = true,
);

sim_force_baryon = Simulation(baryon_particles;
    GravitySolver,
    pids,
);


@info "Initializing MW's halo"
MW_model_halo = Zhao(1.55e7u"Msun/kpc^3", 11.75u"kpc", 1.19, 2.95, 0.95)
MW_ρ_halo_func = (r)->GalacticDynamics.density(MW_model_halo, r)
MW_table_r = collect(0.1:0.1:100) * u"kpc";
MW_ρ_halo = MW_ρ_halo_func.(MW_table_r);
MW_table_mass_halo = 4π * cumul_integrate(MW_table_r, MW_table_r.^2 .* MW_ρ_halo);
MW_table_acc_halo = -C.G * MW_table_mass_halo ./ MW_table_r.^2;
MW_table_pot_halo = -C.G * MW_table_mass_halo ./ MW_table_r;
MW_table_pot_halo[1:end-2] .-= [C.G * 4π * NumericalIntegration.integrate(MW_table_r[i:end], MW_ρ_halo[i:end] .* MW_table_r[i:end]) for i in eachindex(MW_table_r)[1:end-2]];
spl_acc = Spline1D(ustrip.(u"kpc", MW_table_r), ustrip.(u"kpc/Gyr^2", MW_table_acc_halo); k=1)
spl_pot = Spline1D(ustrip.(u"kpc", MW_table_r), ustrip.(u"kpc^2/Gyr^2", MW_table_pot_halo); k=1)


### MW background force without LMC
function MW_bg_force(p::AbstractParticle, t;
    SofteningLength = 0.01u"kpc",
    GravitySolver = Tree(),
    # GravitySolver = DirectSum(),
    sim_force_baryon = sim_force_baryon,
    spl_acc = spl_acc,
)
    acc_DM = spl_acc(ustrip(u"kpc", norm(p.Pos))) * ustrip(normalize(p.Pos)) * u"kpc/Gyr^2"
    acc_b = compute_force(sim_force_baryon, p.Pos, SofteningLength, GravitySolver, CPU())[1]
    return acc = acc_b + acc_DM
end

p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
p1y(sim::Simulation) = find_particle(sim, 1).Pos.y
p1z(sim::Simulation) = find_particle(sim, 1).Pos.z

analysers = Dict(
    "x" => p1x,
    "y" => p1y,
    "z" => p1z,
)



##### LMC
if flag_load_LMC && isfile(file_LMC)
    @info "Loading LMC trajectory from $(file_LMC)"
    df_traj_LMC = DataFrame(CSV.File(file_LMC));
else
    @info "Simulating LMC..."
    ### Vasiliev 2021
    model_LMC = tNFW(5.2629e6u"Msun/kpc^3", 8.5u"kpc" * 1.38^0.6, 85u"kpc" * 1.38^0.6)
    config_LMC = SphericalSystem(STAR, Np_LMC,
        x -> GalacticDynamics.density(model_LMC, x) * 4π * x^2,
    )
    particles_LMC = generate(config_LMC;
        MaxRadius = 85.0u"kpc" * 1.38^0.6 * 1.5,
    );
    particles_LMC.Mass .= 1.38e11u"Msun"/Np_LMC; # make sure the total mass is 1.38e11 Msun


    tidal_initial_pos = PVector(-0.6, -41.3, -27.1, u"kpc") # Galactocentric
    tidal_initial_vel = -PVector(-63.9, -213.8, 206.6, u"km/s") #? Negative for lookback

    test_particles = StructArray([
        Star(uAstro; id = 1),
    ])
    test_particles.Pos[1] = tidal_initial_pos
    test_particles.Vel[1] = tidal_initial_vel

    sim_traj_lookback = Simulation(
        deepcopy(test_particles);
        bgforce = Function[MW_bg_force],
        analysers,
        TimeEnd = 6.0u"Gyr",
        TimeStep = 0.0u"Gyr",  # adaptive if zero.
        # TimeStep = 1e-5u"Gyr",
        TimeBetweenSnapshots = 0.1u"Gyr",
        Realtime = false,
        OutputDir = "output/LMC_traj_lookback_Vasiliev2021",
    );
    run(sim_traj_lookback);


    df_traj_LMC = DataFrame(CSV.File(joinpath("output/LMC_traj_lookback_Vasiliev2021", "analysis.csv")));
end


function MW_bg_force_LMC(p::AbstractParticle, t;
    SofteningLength = 0.01u"kpc",
    GravitySolver = Tree(),
    # GravitySolver = DirectSum(),
    sim_force_baryon = sim_force_baryon,
    spl_acc = spl_acc,
    df_traj = df_traj_LMC,
    particles_LMC = particles_LMC,
)
    acc_DM = spl_acc(ustrip(u"kpc", norm(p.Pos))) * ustrip(normalize(p.Pos)) * u"kpc/Gyr^2"
    acc_b = compute_force(sim_force_baryon, p.Pos, SofteningLength, GravitySolver, CPU())[1]

    if ustrip(u"Gyr", t) <= minimum(df_traj.time)
        id = 1
    else
        id = findfirstvalue(df_traj.time, ustrip(u"Gyr", t))
    end

    particles_LMC_t = deepcopy(particles_LMC)
    particles_LMC_t.Pos .+= PVector(df_traj.x[id], df_traj.y[id], df_traj.z[id]) * u"kpc"
    acc_LMC = AstroNbodySim.compute_force_at_point(p.Pos, particles_LMC_t, C.G, SofteningLength)  #? DirectSum for LMC particles
    return acc = acc_b + acc_DM + acc_LMC
end


@showprogress for i in eachindex(df_MW_satellites.Galaxy)
    @info "Simulating galaxy $(df_MW_satellites.Galaxy[i]) without LMC"
    tidal_initial_pos = PVector(df_MW_satellites.X[i], df_MW_satellites.Y[i], df_MW_satellites.Z[i]) # Galactocentric
    tidal_initial_vel = -PVector(df_MW_satellites.v_X[i], df_MW_satellites.v_Y[i], df_MW_satellites.v_Z[i])  #? Negative for lookback
    @info tidal_initial_pos
    @info tidal_initial_vel

    test_particles = StructArray([
        Star(uAstro; id = 1),
    ])
    test_particles.Pos[1] = tidal_initial_pos
    test_particles.Vel[1] = tidal_initial_vel

    sim_traj_lookback = Simulation(
        deepcopy(test_particles);
        bgforce = Function[MW_bg_force],
        analysers,
        TimeEnd,
        TimeStep,
        TimeBetweenSnapshots,
        Realtime = false,
        OutputDir = joinpath(@__DIR__, "output/traj_lookback_$(df_MW_satellites.Galaxy[i])_no_LMC"),
    );
    run(sim_traj_lookback);


    @info "Simulating galaxy $(df_MW_satellites.Galaxy[i]) with LMC"
    tidal_initial_pos = PVector(df_MW_satellites.X[i], df_MW_satellites.Y[i], df_MW_satellites.Z[i]) # Galactocentric
    tidal_initial_vel = -PVector(df_MW_satellites.v_X[i], df_MW_satellites.v_Y[i], df_MW_satellites.v_Z[i])  #? Negative for lookback

    test_particles = StructArray([
        Star(uAstro; id = 1),
    ])
    test_particles.Pos[1] = tidal_initial_pos
    test_particles.Vel[1] = tidal_initial_vel

    sim_traj_lookback = Simulation(
        deepcopy(test_particles);
        bgforce = Function[MW_bg_force_LMC],
        analysers,
        TimeEnd,
        TimeStep,
        TimeBetweenSnapshots,
        Realtime = false,
        OutputDir = joinpath(@__DIR__, "output/traj_lookback_$(df_MW_satellites.Galaxy[i])_with_LMC"),
    );
    run(sim_traj_lookback);

    #TODO release memory of simulations
end

