# In this example, we generate and simulate a plummer cluster
@info "Loading example 03-plummer"

#=
cd("AstroNbodySim.jl/examples/03-plummer")
=#

# This is a heavy simulation, we would run it in parallel
using AstroNbodySim, PhysicalParticles, Unitful, UnitfulAstro
astro()
mkpathIfNotExist("output")

# AstroIC supports initial condition sampling
using AstroIC
using PhysicalTrees

## First define a config. Keywords are necessary since the config type is immutable
config = PlummerStarCluster(
    collection = STAR,
    NumSamples = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass = 1.0e5u"Msun",
    model = AstroIC.Newton(),
)

## Now generate particles. MaxRadius restricts the sampling region.
particles = generate(config, MaxRadius = 0.1u"kpc");

TimeEnd = 0.1u"Gyr"
TimeStep = 0.000004u"Gyr"
TimeBetweenSnapshots = 0.0005u"Gyr"

# In default, AstroNbodySim runs in master-to-worker mode, 
# which means the master controls the overall planning but not doing the concrete work
@info "Testing master-to-worker mode"
ds = Simulation(
    deepcopy(particles);
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep,
    OutputDir = "output/Plummer-MTW-DirectSum",
);

ds_adapt = Simulation(
    deepcopy(particles);
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep = 0.0u"Gyr",
    OutputDir = "output/Plummer-MTW-DirectSumAdaptive",
);

ts = Simulation(
    deepcopy(particles);
    #octreeconfig = OctreeConfig(length(particles), MaxTopnode = 20000),
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep,
    GravitySolver = Tree(),
    OutputDir = "output/Plummer-MTW-Tree",
);

ts_adapt = Simulation(
    deepcopy(particles);
    #octreeconfig = OctreeConfig(length(particles), MaxTopnode = 20000),
    TimeEnd,
    TimeBetweenSnapshots,
    TimeStep = 0.0u"Gyr",
    GravitySolver = Tree(),
    OutputDir = "output/Plummer-MTW-TreeAdaptive",
);

suggest_softlen!(ds)
suggest_softlen!(ds_adapt)
suggest_softlen!(ts)
suggest_softlen!(ts_adapt)

run(ds)
run(ds_adapt)
run(ts)
run(ts_adapt)

# Plots
using AstroIO

@info "Loading AstroPlot"
using AstroPlot
using AstroPlot.ColorSchemes
using Colors
using GLMakie
using Printf
using FFMPEG

function plotfigs(ds, ds_adapt, ts, ts_adapt, mode::AbstractString)
    @info "Plotting radii"
    ScaleScene, ScaleLayout = layoutscene(resolution = (800, 800))
    LagrangeScene, LagrangeLayout = layoutscene(resolution = (800, 700))

    colors = ColorSchemes.tab10.colors

    AS1 = ScaleLayout[1,1] = GLMakie.Axis(ScaleScene, title = "Direct Sum const dt")
    AS2 = ScaleLayout[1,2] = GLMakie.Axis(ScaleScene, title = "Direct Sum adaptive dt")
    AS3 = ScaleLayout[2,1] = GLMakie.Axis(ScaleScene, title = "Tree const dt")
    AS4 = ScaleLayout[2,2] = GLMakie.Axis(ScaleScene, title = "Tree adaptive dt")

    AL1 = LagrangeLayout[1,1] = GLMakie.Axis(LagrangeScene, title = "Direct Sum const dt")
    AL2 = LagrangeLayout[1,2] = GLMakie.Axis(LagrangeScene, title = "Direct Sum adaptive dt")
    AL3 = LagrangeLayout[2,1] = GLMakie.Axis(LagrangeScene, title = "Tree const dt")
    AL4 = LagrangeLayout[2,2] = GLMakie.Axis(LagrangeScene, title = "Tree adaptive dt")

    plot_radii!(AS1, LagrangeScene, AL1, LagrangeLayout, ds.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(); colors, times = collect(0.0:0.0005:0.1) * u"Gyr", legend=false)
    mv("radii.csv", "output/$(mode)-DirectSum-radii.csv", force = true)

    plot_radii!(AS2, LagrangeScene, AL2, LagrangeLayout, ds_adapt.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(); colors, times = collect(0.0:0.0005:0.1) * u"Gyr", legend=false)
    mv("radii.csv", "output/$(mode)-DirectSumAdaptive-radii.csv", force = true)

    plot_radii!(AS3, LagrangeScene, AL3, LagrangeLayout, ts.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(); colors, times = collect(0.0:0.0005:0.1) * u"Gyr", legend=false)
    mv("radii.csv", "output/$(mode)-Tree-radii.csv", force = true)

    plot_radii!(AS4, LagrangeScene, AL4, LagrangeLayout, ts_adapt.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(); colors, times = collect(0.0:0.0005:0.1) * u"Gyr", legend=false)
    mv("radii.csv", "output/$(mode)-TreeAdaptive-radii.csv", force = true)

    #LagrangeLayout[:,3] = GLMakie.Legend(
    #    LagrangeScene,
    #    LagrangeScene.children[1].plots[2:end],
    #    ["10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%"];
    #    tellheight = false,
    #    tellwidth = false,
    #    halign = :left,
    #    valign = :center,
    #    margin = (0,0,0,0),
    #)

    colsize!(LagrangeLayout, 1, Relative(0.47))
    colsize!(LagrangeLayout, 2, Relative(0.47))
    rowsize!(LagrangeLayout, 1, Relative(0.5))
    rowsize!(LagrangeLayout, 2, Relative(0.5))

    supertitle = ScaleLayout[0,:] = Label(ScaleScene, "Scale Radius")
    Makie.save("output/$(mode)-ScaleRadius.png", ScaleScene)

    supertitle = LagrangeLayout[0,:] = Label(LagrangeScene, "Lagrange Radii")
    Makie.save("output/$(mode)-LagrangianRadii.png", LagrangeScene)
    

    #@info "Plotting positions"
    #plot_positionslice(ds.config.output.dir,       "snapshot_", collect(0:200), ".gadget2", gadget2(), dpi = 300, resolution = (800,800),
    #                   xlims = (-0.05, +0.05), ylims = (-0.05, +0.05), times = collect(0.0:0.0005:0.1) * u"Gyr")
    #plot_positionslice(ds_adapt.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(), dpi = 300, resolution = (800,800),
    #                   xlims = (-0.05, +0.05), ylims = (-0.05, +0.05), times = collect(0.0:0.0005:0.1) * u"Gyr")
    #plot_positionslice(ts.config.output.dir,       "snapshot_", collect(0:200), ".gadget2", gadget2(), dpi = 300, resolution = (800,800),
    #                   xlims = (-0.05, +0.05), ylims = (-0.05, +0.05), times = collect(0.0:0.0005:0.1) * u"Gyr")
    #plot_positionslice(ts_adapt.config.output.dir, "snapshot_", collect(0:200), ".gadget2", gadget2(), dpi = 300, resolution = (800,800),
    #                   xlims = (-0.05, +0.05), ylims = (-0.05, +0.05), times = collect(0.0:0.0005:0.1) * u"Gyr")
   

    ## Animation
    #fps = 25
    #loop = 0
    #ffmpeg_exe(`-v error -framerate $fps -loop $loop -i $(ds.config.output.dir      )/pos_%04d.png output/$(mode)-DirectSum.mp4`)
    #ffmpeg_exe(`-v error -framerate $fps -loop $loop -i $(ds_adapt.config.output.dir)/pos_%04d.png output/$(mode)-DirectSumAdaptive.mp4`)
    #ffmpeg_exe(`-v error -framerate $fps -loop $loop -i $(ts.config.output.dir      )/pos_%04d.png output/$(mode)-Tree.mp4`)
    #ffmpeg_exe(`-v error -framerate $fps -loop $loop -i $(ts_adapt.config.output.dir)/pos_%04d.png output/$(mode)-TreeAdaptive.mp4`)
end

plotfigs(ds, ds_adapt, ts, ts_adapt, "MTW")


# write success flag for shell
success = open("output/success", "w")
close(success)

# Clear memory
AstroNbodySim.clear()