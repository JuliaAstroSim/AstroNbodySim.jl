using Pkg

@info "Please follow the development guides in README before running examples"
@info "Packages in JuliaAstroSim have to be installed manually. It is recommanded to use the latest master branch"

Pkg.add("Unitful")
Pkg.add("UnitfulAstro")

Pkg.add("Colors")
Pkg.add("ColorSchemes")
Pkg.add("FFMPEG")
Pkg.add("GLMakie")
Pkg.add("UnicodePlots")

Pkg.add("CUDA")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("FileIO")
Pkg.add("AstroLib")
Pkg.add("StructArrays")
Pkg.add("StaticArrays")
Pkg.add("Measurements")
Pkg.add("Zygote")