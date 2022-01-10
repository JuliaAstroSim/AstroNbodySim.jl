"""
Compile with:
julia --project=docs/ --color=yes docs/make.jl

Generate key:
DocumenterTools.genkeys(user="JuliaAstroSim", repo="git@github.com:JuliaAstroSim/AstroNbodySim.jl.git")
"""

using Documenter

using AstroNbodySim

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

makedocs(
    modules = [AstroNbodySim],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaastrosims.github.io/AstroNbodySim.jl/dev/",
        assets = ["assets/alpha_small.ico"],
        analytics = "UA-153693590-1",
        highlights = ["llvm", "yaml"],
    ),
    clean = false,
    sitename = "AstroNbodySim.jl",
    authors = "islent",
    linkcheck = false,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial/basic.md",
        ],
        "Manual" => Any[
            "manual/DirectSum.md",
        ],
        "Examples" => Any[
            "examples/01-Binary.md",
            "examples/02-autodiff-bg.md",
            "examples/03-Measurements.md",
            "examples/04-Bigfloat.md",
            "examples/05-plummer.md",
            "examples/06-GalaxyCollision.md",
            "examples/07-TDEcluster.md",
            "examples/08-SolarSystem.md",
        ],
        "Library" => Any[
            "lib/Types.md",
            "lib/Methods.md",
        ],
        "CITATION" => "citation.md",
        "Contributing" => "contributing.md",
        "中文教程" => Any[
            "ChineseTutorial/basic.md",
        ],
        "中文手册" => Any[
            "ChineseManual/DirectSum.md",
        ],
    ],
    #workdir = joinpath(@__DIR__, "../test"),
    #strict = !("strict=false" in ARGS),
    #doctest = ("doctest=only" in ARGS) ? :only : true,
)

deploydocs(
    repo = "github.com/JuliaAstroSim/AstroNbodySim.jl.git",
    target = "build",
    devbranch = "main",
)