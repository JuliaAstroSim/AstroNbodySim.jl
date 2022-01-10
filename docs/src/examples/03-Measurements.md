# 03 Measurements: Uncertainty propagation

In this example, we demonstrate how to:
- Compute with `Measurement`
- Propagate uncertainties from initial conditions
- Plot orbit uncertainties

```julia
# The basic setup is same with example-01
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

p1x(sim::Simulation) = find_particle(sim, 1).Pos.x
p1y(sim::Simulation) = find_particle(sim, 1).Pos.y

analysers = Dict(
    "x" => p1x,
    "y" => p1y,
)


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
```

![Uncertainty of elliptic orbit](https://github.com/JuliaAstroSim/AstroNbodySim.jl/tree/main/docs/src/examples/pics/examples/01-binary/Uncertainty%20of%20elliptic%20orbit.png "Uncertainty of elliptic orbit")
