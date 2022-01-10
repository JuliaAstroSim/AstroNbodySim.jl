# 01 Binary orbit

In this example, we demonstrate how to:
- Output simulation information on the run
- Estimate properties of binary orbit and config the simulation
- Manually generate binary stars
- Plot orbits with gapping
- Add a background force field

## Circular binary

```julia
using AstroIO
using AstroPlot
using AstroPlot.ColorSchemes
using Colors
using CSV, DataFrames

using Plots
pyplot() # Need PyPlot to be installed

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

    p = Plots.plot(plot_x, plot_y,
        aspect_ratio = 1, legend = nothing, xlabel = "x [kpc]", ylabel = "y [kpc]",
        size = resolution,
    )
    savefig(p, "output/" * title * ".png")
    @info "Total orbits: $(div(flip, 2))"
    @info "Steps per orbit: $(length(y) / div(flip, 2))"
    @info "Figure has been saved to " * title * ".png"
    return df
end
```

`analysers::Dict` has full access to simulation data.
The key values are functions that return one printable value, and they will be written in `analysis.csv`.
Here we find one of the particle and output its x- and y- coordinates.

```julia
p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
p1y(sim::Simulation) = find_particle(sim, 1).Pos.y

analysers = Dict(
    "x" => p1x,
    "y" => p1y,
)
```

Now we estimate the orbit properties and set up the simulation

```julia
G = Constant(uAstro).G
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
```

## Massless circular binary with background force field

```julia
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
```