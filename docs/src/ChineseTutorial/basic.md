# 基础使用

在着手模拟之前，建议首先熟悉 `AstroNbodySim.jl` 的依赖包

整个项目有较为完善的注释，在开发和使用过程中如果不熟悉接口的参数和关键字，可以在 `REPL` 中使用 `help?>` 快速获得帮助，比如：
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

## 带单位的矢量和粒子类型

```@repl basic
# PhysicalParticles 定义了矢量和粒子类型
using PhysicalParticles

# 将默认（换算）单位设置为天体物理单位
astro()

# using UnitfulAstro 十分必要，如果在使用天体物理单位的话
using UnitfulAstro

# 定义带单位的矢量
a = PVector(3.0u"kpc", 4.0u"kpc", 12.0u"kpc")
b = PVector(1.0, 1.0, 1.0, u"kpc")

c = PVector2D()
d = PVector2D(0.0, 1.0)

# 基本矢量运算
a * b
c + d

norm(a)
normalize(a)

# 矢量的数组
points = rand(PVector{Float64}, 5) * u"kpc"
p = randn_pvector2d(5)

mean(p)
PhysicalParticles.center(p)
median(p)

# 定义粒子
particles = [Star(uAstro, id = i) for i in 1:5]

# 逐一修改粒子坐标
assign_particles(particles, :Pos, points)

# 平均坐标
average(particles, :Pos)

# 质心坐标
assign_particles(particles, :Mass, rand(5) * u"Msun")
averagebymass(particles, :Pos)

# StructArray 可以更加高效地修改粒子数据
StructArray(particles)

# 可以直接创建 StructArray
s = StructArray(Star(uAstro, id = i) for i in 1:5)

# StructArray 同样支持其他函数，比如 assign_particles, averagebymass 等
# It is much more convenient to use dot operations
s.Pos .= points
s.Mass .= rand(5) * u"Msun"
s
```

## 生成初始条件

```@repl basic
using AstroIC
config = PlummerStarCluster()
particles = generate(config)
```

## 可视化

```@example basic
using AstroPlot
scene, layout = plot_makie(particles)
scene
```

## Snapshot 文件读写

```@repl basic
using AstroIO
if !isdir("output/")
    mkpath("output/")
end

write_csv("output/basic.csv", particles)
write_jld("output/basic.jld2", particles)
write_gadget2("output/basic.gadget2", particles) # This would generate a header automatically

# 可以从粒子数据生成 Gadget2 header
header = HeaderGadget2(particles)
header.time = 0.1   # Gyr

write_gadget2("output/basicwithheader.gadget2", header, particles, uGadget2) # write in Gadget2 units (default)

# 读取刚刚写入的文件
h, d = read_gadget2("output/basic.gadget2", uAstro)
d = read_jld("output/basic.jld2")
```