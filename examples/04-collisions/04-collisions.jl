# In this example, we compare the galaxy collision example with Gadget-2
# and collide galaxies to explore a parameter space

#=
cd("AstroNbodySim.jl/examples/04-collisions")
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

header, d = read_gadget2(joinpath(@__DIR__, "galaxy_littleendian.dat"), uAstro, type=Star);  # 20000 disk and 40000 halo particles
#d64 = convert.(Float64, d)
#data = convert.(Int, d64)

#=
# Gadget run

rm("output/gadget", force = true, recursive = true)
mkdir("output/gadget")
cd("output")

# compile
cp("../../Gadget2", "Gadget2", force = true)
cd("Gadget2")
write_gadget2_makefile("Makefile")

# Make sure that gravity constant is the same
G_cgs = ustrip(uconvert(u"cm^3/g/s^2", Constant().G))
alter_line("allvars.h", "#define  GRAVITY",
"#define  GRAVITY           $(G_cgs)   /*!< Gravitational constant (in cgs units) */")

result = run(`make -j 4`)
cd("../")

result = run(`mpirun -np 4 Gadget2/Gadget2 ../galaxy.param`)
cd("../")

gpu = Simulation(
    deepcopy(data);
    units = uAstro,
    TimeEnd,
    TimeBetweenSnapshots,
    ErrTolTimestep = 0.025,
    OutputDir = "output/Collision-DirectSumAdaptiveGPU",
    Realtime = false,
    device = GPU(),
);
set_softlen!(gpu, 0.4u"kpc")
=#

TimeEnd = 3.0f0u"Gyr"
TimeBetweenSnapshots = 0.03f0u"Gyr"

gpu = Simulation(
    deepcopy(d);
    floattype = Float32,
    units = uAstro,
    TimeEnd,
    TimeBetweenSnapshots,
    OutputDir = "output/Collision-DirectSumAdaptiveGPU",
    constants = Constant(Float32, uAstro),
    ZeroValues = ZeroValue(Float32, uAstro),
    ForceSofteningTable = [0.4f0u"kpc" for i in 1:6],
    Realtime = false,
    device = GPU(),
);
run(gpu)


# Plot
plot_positionslice(gpu.config.output.dir, "snapshot_", collect(0:100), ".gadget2", gadget2(),
    size = (400,400),
    xlims = (-200.0, +200.0), ylims = (-200.0, +200.0),
    times = collect(0.0:0.03:3.0) * u"Gyr",
    collection = DISK,
    markersize = 2.0,
)
plt = mosaic(gpu.config.output.dir, "pos_", collect(1:9:100), ".png"; fillvalue = 0.5, npad = 3, ncol = 4, rowmajor = true);
save("output/mosaic-collision-DirectSumAdaptiveGPU.png", plt)

# write success flag for shell
success = open("output/success", "w")
close(success)

# Clear memory
AstroNbodySim.clear()