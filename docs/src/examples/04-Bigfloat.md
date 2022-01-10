# 04 BigFloat

In this example, we demonstrate how to compute with `BigFloat`.

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
```