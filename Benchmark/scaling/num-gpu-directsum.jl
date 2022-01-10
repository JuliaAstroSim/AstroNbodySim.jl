@info "Loading"

using AstroNbodySim
using AstroIC
using PhysicalParticles, UnitfulAstro

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
    device = GPU(),
    OutputDir = "output/GPU-DirectSum-$(length(data))",
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
        device = GPU(),
        OutputDir = "output/GPU-DirectSum-$(length(d))",
    )
    run(ds)
end