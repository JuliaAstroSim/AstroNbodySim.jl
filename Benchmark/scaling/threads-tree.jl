@info "Loading"
println("Number of threads: ", Threads.nthreads())

using AstroNbodySim
using AstroIO
using PhysicalParticles, UnitfulAstro

TimeEnd = 0.00005u"Gyr"
TimeStep = 0.00001u"Gyr"
TimeBetweenSnapshots = TimeEnd

header, data = read_gadget2("output/scaling.gadget2", uAstro, type = Star);
d64 = convert.(Float64, data)
d = convert.(Int64, d64)
ts = Simulation(
    d;
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep,
    GravitySolver = Tree(),
    OutputDir = "output/threads-Tree-$(Threads.nthreads())",
)
run(ts)