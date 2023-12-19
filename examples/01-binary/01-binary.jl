@info "Loading example 01-binary"

#=
cd("AstroNbodySim.jl/examples/01-binary")
include("01-binary.jl")
=#

using AstroIO
using AstroPlot
using AstroPlot.ColorSchemes
using Colors
using CSV, DataFrames

using Plots
# pyplot() # Need PyPlot to be installed

using AstroNbodySim
using PhysicalParticles
using Unitful, UnitfulAstro
astro()
mkpathIfNotExist("output")

function plot_orbit_with_gap(sim::Simulation, title::String, gap::Int = 1; resolution = (800,800), kw...)
    df = DataFrame(CSV.File(joinpath(sim.config.output.dir, "analysis.csv")))
    x = df.x
    y = df.y

    title *= " of every $(gap) orbit(s)"
    
    count = 0
    flip = 0
    Head = 1
    last_sign = sign(y[1])

    plot_x = []
    plot_y = []
    for i in 1:length(y)
        @inbounds new_sign = sign(y[i])
        if new_sign != last_sign    # half cycle
            flip += 1
            last_sign = new_sign
            
            if flip % 2 == 0        # full cycle
                count += 1
                if count % gap == 0 # plot this cycle
                    push!(plot_x, x[Head:i])
                    push!(plot_y, y[Head:i])
                    count = 0       # Reset gap
                end
                Head = i            # start of cycle
            end
        end
        # do nothing
    end

    # We have to use Plots to avoid autolimitaspect StackOverflowError:
    # https://github.com/JuliaPlots/Makie.jl/issues/1525
    p = Plots.plot(plot_x, plot_y,
        #title = title,
        aspect_ratio = 1, legend = nothing, xlabel = "x [kpc]", ylabel = "y [kpc]",
        size = resolution,
    )
    savefig(p, "output/" * title * ".png")
    @info "Total orbits: $(div(flip, 2))"
    @info "Steps per orbit: $(length(y) / div(flip, 2))"
    @info "Figure has been saved to " * title * ".png"
    return df
end

# 
p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
p1y(sim::Simulation) = find_particle(sim, 1).Pos.y

analysers = Dict(
    "x" => p1x,
    "y" => p1y,
)

G = Constant(uAstro).G


### Binary Circular
mass = 1.0e10u"Msun"
D = 2.0u"kpc"

binary_orbit_period(mass, D) = 2 * pi * sqrt(0.5 * D^3 / G / mass)

binary_rot_vel(mass, D) = sqrt(G * mass / D / 2.0)

vel = binary_rot_vel(mass, D)
T = binary_orbit_period(mass, D)

@info "Estimated Period: $(T)"
@info "Estimated Rotation Vel: $(vel)"

# Define the two particles
data = StructArray(Star(uAstro, id = i) for i in 1:2);
data.Pos[1] = PVector(+0.5D, 0.0u"kpc", 0.0u"kpc");
data.Pos[2] = PVector(-0.5D, 0.0u"kpc", 0.0u"kpc");
data.Vel[1] = PVector(0.0u"kpc/Gyr", +vel, 0.0u"kpc/Gyr");
data.Vel[2] = PVector(0.0u"kpc/Gyr", -vel, 0.0u"kpc/Gyr");
data.Mass .= mass;

# Define simulations
function compute_dt(a::Number; SofteningLength = 0.0001u"kpc", ErrTolTimestep = 0.025)
    return sqrt(2 * ErrTolTimestep * SofteningLength / abs(a))
end

TimeStep = compute_dt(G*mass/D^2)
TimeEnd = T * 1.01
TimeBetweenSnapshots = TimeEnd

circular = Simulation(
    deepcopy(data);
    analysers,
    TimeEnd,
    TimeStep,
    TimeBetweenSnapshots,
    OutputDir = "output/BinaryCircular",
);
run(circular)
plot_orbit_with_gap(circular, "Binary Circular")




### Now we add a background force field and simulate two massless particles
mass = 0.5e10u"Msun"
R = 1.0u"kpc"

binary_orbit_period(mass, R) = 2 * pi * sqrt(R^3 / G / mass)
binary_rot_vel(mass, R) = sqrt(G * mass / R)

vel = binary_rot_vel(mass, R)
T = binary_orbit_period(mass, R)

@info "Estimated Period: $(T)"
@info "Estimated Rotation Vel: $(vel)"

# The only difference is that the particles have no mass, so they cannot produce any mutual forces
data = StructArray(Star(uAstro, id = i) for i in 1:2);
data.Pos[1] = PVector(+R, 0.0u"kpc", 0.0u"kpc");
data.Pos[2] = PVector(-R, 0.0u"kpc", 0.0u"kpc");
data.Vel[1] = PVector(0.0u"kpc/Gyr", +vel, 0.0u"kpc/Gyr");
data.Vel[2] = PVector(0.0u"kpc/Gyr", -vel, 0.0u"kpc/Gyr");

# Define the background force field
attractor(p::AbstractParticle) = -G * mass / (p.Pos * p.Pos) * normalize(p.Pos) / 1.0u"kpc"
bgforce = Function[attractor]

TimeEnd = T * 1.01
TimeStep = compute_dt(G*mass/R^2)
TimeBetweenSnapshots = TimeEnd

bg = Simulation(
    deepcopy(data);
    bgforce,
    analysers,
    TimeEnd,
    TimeStep,
    TimeBetweenSnapshots,
    OutputDir = "output/bgforce",
);
run(bg)
plot_orbit_with_gap(bg, "Binary Circular with bgforce")




### Elliptic orbit precision
R = 1.0u"kpc"
mass = 1.0e8u"Msun"

e = 0.9

vel = sqrt((1.0-e)*G*mass/R)
@info "Velocity at ap: $(vel)"

a = ellipticSemiMajor(G, mass, R, vel)
@info "Length of semi-major axis: $(a)"

T = ellipticPeriod(G, mass, a)
@info "Period of elliptical orbit: $(T)"

data = StructArray(Star(uAstro, id = i) for i in 1:2)
data.Mass[2] = mass

data.Pos[1] = PVector(R, 0.0u"kpc", 0.0u"kpc")
data.Vel[1] = PVector(0.0u"kpc/Gyr", vel, 0.0u"kpc/Gyr")


NumOrbits = 200
TimeEnd = NumOrbits * T * 1.001
TimeBetweenSnapshots = TimeEnd

# Adaptive timesteps
elliptic_adapt = Simulation(
    deepcopy(data);
    analysers,
    TimeEnd,
    ErrTolTimestep = 0.05,
    TimeBetweenSnapshots,
    ForceSofteningTable = [0.001u"kpc" for i in 1:6],
    OutputDir = "output/EllipticOrbitAdapt",
)
run(elliptic_adapt)
df = plot_orbit_with_gap(elliptic_adapt, "Elliptic Orbit (adaptive)", 40, resolution = (600,300))

# Fixed timesteps
Δt = df.time[2:end-1] .- df.time[1:end-2]
elliptic_const = Simulation(
    deepcopy(data);
    analysers,
    TimeEnd,
    TimeStep = minimum(Δt) * u"Gyr",
    ErrTolTimestep = 0.1,
    TimeBetweenSnapshots,
    ForceSofteningTable = [0.01u"kpc" for i in 1:6],
    OutputDir = "output/EllipticOrbitConst",
)
run(elliptic_const)
df = plot_orbit_with_gap(elliptic_const, "Elliptic Orbit (const)", 40, resolution = (600,300))



### Measurements
# We use the same elliptic orbit
#TODO: plot x-axis and y-axis uncertainties at the same time
using Measurements
G = Constant(Measurement, uAstro).G

mass = measurement(1.0e8u"Msun", 1.0e6u"Msun")
R = measurement(1.0u"kpc", 0.01u"kpc")

e = 0.9

vel = sqrt((1.0-e)*G*mass/R)
@info "Velocity at ap: $(vel)"

a = ellipticSemiMajor(G, mass, R, vel)
@info "Length of semi-major axis: $(a)"

T = ellipticPeriod(G, mass, a)
@info "Period of elliptical orbit: $(T)"

# Define the two particles
data = StructArray(Star(Measurement, Int, uAstro, id = i) for i in 1:2)
data.Pos[1] = PVector(R, measurement(0.0u"kpc"), measurement(0.0u"kpc"))
data.Vel[1] = PVector(measurement(0.0u"kpc/Gyr"), +vel, measurement(0.0u"kpc/Gyr"))

data.Mass[2] = mass


TimeEnd = Measurements.value(T) * 0.95
TimeBetweenSnapshots = TimeEnd

m = Simulation(
    deepcopy(data);
    analysers,
    TimeEnd,
    TimeBetweenSnapshots,
    OutputDir = "output/Measurements",
    constants = Constant(Measurement, uAstro),
    ZeroValues = ZeroValue(Measurement, uAstro),
    ForceSofteningTable = measurement.([1.0e-4u"kpc" for i in 1:6]),
)
run(m)

# plot the orbit with uncertainties
function plot_orbit_with_uncertainties(sim::Simulation, title::String;
    resolution = (800,450),
    xlabel = "x [kpc]", ylabel = "y [kpc]",
)
    df = DataFrame(CSV.File(joinpath(sim.config.output.dir, "analysis.csv")))
    xm = measurement.(df.x)
    ym = measurement.(df.y)

    x = Measurements.value.(xm)
    y = Measurements.value.(ym)

    # find the flip one, to avoid wrong color of covered error bars
    s = sign(y[1])
    flip = 0
    for i in eachindex(y)
        if s != sign(y[i])
            flip = i
            break
        end
    end

    xerror = Measurements.uncertainty.(xm)
    yerror = Measurements.uncertainty.(ym)

    p = Plots.plot(x[1:flip], y[1:flip];
        ribbon = yerror[1:flip], legend=nothing,
        xlabel, ylabel, title,
        size = resolution,
        aspect_ratio = 1,
    )
    Plots.plot!(p, x[flip+1:end], y[flip+1:end]; ribbon = yerror[flip+1:end])
    Plots.plot!(p, [1.0; 0.92], [0.1,0.15],arrow=(2,1.0))
    savefig(p, "output/" * title * ".png")
    return DataFrame(time = df.time, x = xm, y = ym)
end
plot_orbit_with_uncertainties(m, "Uncertainty of elliptic orbit", resolution = (800,450))

# autodiff?
Measurements.derivative(m.simdata.Pos[1].y.val, m.simdata.Mass[2].val)

Measurements.uncertainty_components(m.simdata.Pos[1].y.val)



### Bigfloat
# Here we only demonstrate how to construct a simulation using BigFloat
astro()
G = Constant(BigFloat, uAstro).G

mass = BigFloat("1.0e10")u"Msun"
D = BigFloat("2.0")u"kpc"

binary_orbit_period(mass, D) = 2 * pi * sqrt(0.5 * D^3 / G / mass)

binary_rot_vel(mass, D) = sqrt(G * mass / D / 2.0)

vel = binary_rot_vel(mass, D)
T = binary_orbit_period(mass, D)

@info "Estimated Period: $(T)"
@info "Estimated Rotation Vel: $(vel)"

data = StructArray(Star(BigFloat, Int, uAstro, id = i) for i in 1:2)
data.Pos[1] = PVector(+0.5*BigFloat(D), BigFloat(0.0)*u"kpc", BigFloat(0.0)*u"kpc")
data.Pos[2] = PVector(-0.5*BigFloat(D), BigFloat(0.0)*u"kpc", BigFloat(0.0)*u"kpc")
data.Vel[1] = PVector(BigFloat(0.0)*u"kpc/Gyr", +BigFloat(vel), BigFloat(0.0)*u"kpc/Gyr")
data.Vel[2] = PVector(BigFloat(0.0)*u"kpc/Gyr", -BigFloat(vel), BigFloat(0.0)*u"kpc/Gyr")
data.Mass .= mass

TimeEnd = 0.030u"Gyr"
TimeStep = 0.00002u"Gyr"
TimeBetweenSnapshots = 0.0005u"Gyr"

big = Simulation(
    deepcopy(data);
    TimeEnd,
    TimeBetweenSnapshots, 
    TimeStep,
    OutputDir = "output/BigFloat",
    constants = Constant(BigFloat, uAstro),
    ZeroValues = ZeroValue(BigFloat, uAstro),
    ForceSofteningTable = [BigFloat("1.0e-4")*u"kpc" for i in 1:6],
)
compute_force(big)
compute_potential(big)

display(big.simdata.Acc)
display(big.simdata.Potential)

run(big)
display(big.simdata.Pos)



# write success flag for shell
success = open("output/success", "w")
close(success)

# Clear memory
AstroNbodySim.clear()