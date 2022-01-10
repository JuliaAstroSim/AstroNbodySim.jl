@info "Loading"
println("Number of threads: ", Threads.nthreads())

using AstroNbodySim
using PhysicalParticles, UnitfulAstro
using AstroIC

TimeEnd = 0.00005u"Gyr"
TimeStep = 0.00001u"Gyr"
TimeBetweenSnapshots = TimeEnd

@info "warm up"
data = generate(PlummerStarCluster(NumSamples = 100))
warm = Simulation(
    data;
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep,
    OutputDir = "output/threads-DirectSum-$(length(data))",
)
run(warm)

NumData = parse.(Int, ARGS)
for N in NumData
    d = generate(PlummerStarCluster(NumSamples = N))
    ds = Simulation(
        d;
        TimeEnd,
        TimeBetweenSnapshots,
        TimeStep,
        OutputDir = "output/threads-DirectSum-$(length(d))",
    )
    run(ds)
end