##### ParameterSpace test #####

#=
cd("AstroNbodySim/examples/04-collisions")
include("04-collisions.jl")
=#


using Unitful, UnitfulAstro
using Distributed
using PhysicalParticles
using AstroIC, AstroIO
using AstroNbodySim
using AstroPlot
using GLMakie
using FFMPEG
using Printf

mkpathIfNotExist("output")

#TODO Change to disk galaxies
config = PlummerStarCluster(NumSamples = 1000)
galaxy1 = generate(config);
galaxy2 = generate(config);
galaxy3 = generate(config);
galaxy4 = generate(config);

setpos(galaxy1, PVector(-0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy2, PVector(+0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy3, PVector(0.0, +0.1, 0.0, u"kpc"))
setpos(galaxy4, PVector(0.0, -0.1, 0.0, u"kpc"))

setvel(galaxy1, PVector(+0.4, +0.8, 0.0, u"kpc/Gyr"))
setvel(galaxy2, PVector(-0.4, -1.0, -0.2, u"kpc/Gyr"))
setvel(galaxy3, PVector(+1.0, -0.4, 0.4, u"kpc/Gyr"))
setvel(galaxy4, PVector(-1.0, +0.4, -0.2, u"kpc/Gyr"))

data = deepcopy(galaxy1);
append!(data, galaxy2);
append!(data, galaxy3);
append!(data, galaxy4);
data.ID .= 1:4000;


TimeEnd = 1.0u"Gyr"
TimeBetweenSnapshots = 0.002u"Gyr"

SoftLen = suggest_softlen(galaxy1) / 10
ForceSofteningTable = [SoftLen for i in 1:6]

sim = DirectSumGPUSimulation(
    data,
    config = SimConfig(;
        TimeEnd,
        TimeBetweenSnapshots,
        ForceSofteningTable,
        ErrTolTimestep = 0.01,
        OutputDir = "output/test",
        SaveRestart = false,
    )
);

#=
run(sim)
=#


#=
plot_positionslice("output/test", "snapshot_", collect(0:500), ".gadget2", gadget2(), resolution = (800, 800),
                   xlims = (-0.15,0.15), ylims = (-0.15,0.15), markersize = 1.0)
ffmpeg_exe(`-v error -framerate 25 -loop 0 -i output/test/pos_%04d.png output/test.mp4`)
=#

### Tuning on impact parameters
using ParameterSpace
params = [
    Parameter("angle", 1, pi/12:pi/12:pi/4), # rad
    Parameter("vel", 2, [0.5, 1.0]),       # kpc/Gyr
    #Parameter("mass1"),
    #Parameter("mass2"),
]

function galaxy_colliders(angle::Float64, vel::Number;
        TimeEnd = 0.5u"Gyr",
        TimeBetweenSnapshots = 0.001u"Gyr",
    )
    config = PlummerStarCluster(NumSamples = 2500)
    galaxy1 = generate(config)
    galaxy2 = generate(config)

    vx = vel * cos(angle)
    vy = vel * sin(angle)

    setpos(galaxy1, PVector(-0.1, 0.0, 0.0, u"kpc"))
    setpos(galaxy2, PVector(+0.1, 0.0, 0.0, u"kpc"))

    setvel(galaxy1, PVector(+vx, +vy, 0.0, u"kpc/Gyr"))
    setvel(galaxy2, PVector(-vx, -vy, 0.0, u"kpc/Gyr"))

    data = deepcopy(galaxy1)
    append!(data, galaxy2)

    simlabel = @sprintf("%.2f - %.2f", angle, vel)
    sim = DirectSumGPUSimulation(
        data,
        config = SimConfig(;
            TimeEnd,
            TimeBetweenSnapshots,
            OutputDir = "output/$simlabel",
            SaveRestart = false,
        ),
        Realtime = false,
    )
    run(sim)

    plot_positionslice(sim.config.output.dir, "snapshot_", collect(0:Int(TimeEnd/TimeBetweenSnapshots)), ".gadget2", gadget2(), resolution = (960, 960),
        xlims = (-0.15,0.15), ylims = (-0.15,0.15), markersize = 0.1, times = collect(0.0u"Gyr":TimeBetweenSnapshots:TimeEnd))
    
    ffmpeg_exe(`-v error -framerate 25 -loop 0 -i $(sim.config.output.dir)/pos_%04d.png output/'video - '$(simlabel).mp4`)

    fig = plot_rotationcurve(sim.simdata, savefolder = sim.config.output.dir)
    Makie.save("output/RotationCurve - $simlabel.png", fig)

    Counts = collect(0:Int(TimeEnd/TimeBetweenSnapshots/50):Int(TimeEnd/TimeBetweenSnapshots))[1:end-1]
    fig = mosaicview(sim.config.output.dir, "pos_", Counts, ".png"; fillvalue=0.5, npad=5, ncol=10, rowmajor=true)
    Makie.save("output/mosaic - $simlabel.png", fig)

    # release memory of former simulations
    pop!(AstroNbodySim.registry, sim.id)
    return simlabel # Test
end

#=
df = analyse_function(galaxy_colliders, params)
=#