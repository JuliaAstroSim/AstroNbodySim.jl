# Basic

Before working on a simulation, it is recommanded to be familiar with dependencies of `AstroNbodySim.jl`

The whole project is well documented.
Use `help?` in `REPL` is the fastest way to get help if you forget the arguments or keywords of interfaces.
For example,
```julia
julia> using AstroIO

help?> read_gadget2
search: read_gadget2 read_gadget2_pos read_gadget2_jld

  read_gadget2(filename::AbstractString, units, fileunits = uGadget2; kw...)

  Return a Tuple of header and particle data in snapshot file.
  units is supported by PhysicalParticles: uSI, uCGS, uAstro,
  uGadget2, nothing. fileunits is the internal units in the file,
  and will be converted to units while reading the file.

  Keywords
  ≡≡≡≡≡≡≡≡≡≡

    •  acc::Bool = false : read acceleration data if exist

    •  pot::Bool = false : read potential data if exist
```

## Physical vectors and particles

```@repl basic
# PhysicalParticles defines the basic physical vector and particle type
using PhysicalParticles

# Set the default units to astrophysics. See `Unitful.preferunits` for more info
astro()

# using UnitfulAstro is necessary, if you are using astrophysical units
using UnitfulAstro

# Define a physical vector
a = PVector(3.0u"kpc", 4.0u"kpc", 12.0u"kpc")
b = PVector(1.0, 1.0, 1.0, u"kpc")

c = PVector2D()
d = PVector2D(0.0, 1.0)

# Basic algebra
a * b
c + d

norm(a)
normalize(a)

# generate an array of physical vector
points = rand(PVector{Float64}, 5) * u"kpc"
p = randn_pvector2d(5)

mean(p)
PhysicalParticles.center(p)
median(p)

# Define particles
particles = [Star(uAstro, id = i) for i in 1:5]

# Assign particle positions with array of points
assign_particles(particles, :Pos, points)

# average position
average(particles, :Pos)

# average position by mass
assign_particles(particles, :Mass, rand(5) * u"Msun")
averagebymass(particles, :Pos)

# StructArray is more efficient to manipulate field of struct in an array
StructArray(particles)

# Or construct directly
s = StructArray(Star(uAstro, id = i) for i in 1:5)

# StructArray is still supported by assign_particles, averagebymass, etc.
# It is much more convenient to use dot operations
s.Pos .= points
s.Mass .= rand(5) * u"Msun"
s
```

## Generate initial conditions

```@repl basic
using AstroIC
config = PlummerStarCluster()
particles = generate(config)
```

## Visualization

```@example basic
using AstroPlot
fig = plot_makie(particles)
fig
```

## Snapshot File I/O

```@repl basic
using AstroIO
if !isdir("output/")
    mkpath("output/")
end

write_csv("output/basic.csv", particles)
write_jld("output/basic.jld2", particles)
write_gadget2("output/basic.gadget2", particles) # This would generate a header automatically

# generate a header and modify it
header = HeaderGadget2(particles)
header.time = 0.1   # Gyr

write_gadget2("output/basicwithheader.gadget2", header, particles, uGadget2) # write in Gadget2 units (default)

# now load the files we had just written
h, d = read_gadget2("output/basic.gadget2", uAstro)
d = read_jld("output/basic.jld2")
```