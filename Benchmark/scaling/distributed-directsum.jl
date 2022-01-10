@info "Loading"
using Distributed
println("Number of processes: ", nprocs())

@everywhere using AstroNbodySim
@everywhere using AstroIO
@everywhere using PhysicalParticles, UnitfulAstro

TimeEnd = 0.00005u"Gyr"
TimeStep = 0.00001u"Gyr"
TimeBetweenSnapshots = TimeEnd

header, data = read_gadget2("output/scaling.gadget2", uAstro, type = Star);
d64 = convert.(Float64, data)
d = convert.(Int64, d64)
ds = Simulation(
    d;
    pids = procs(),
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep,
    OutputDir = "output/distributed-DirectSum-$(nprocs())",
)
run(ds)