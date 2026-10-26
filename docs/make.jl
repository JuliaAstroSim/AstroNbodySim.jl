"""
Init with:
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

Compile with:
julia --project=docs/ --color=yes docs/make.jl

Deployment is handled out-of-band by the `julia-actions/julia-docdeploy`
GitHub Action (see `.github/workflows/ManuallyBuildDoc.yml`).
"""

using Documenter

using AstroNbodySim

# The DOCSARGS environment variable can be used to pass additional
# arguments to make.jl.  Useful for CI when tweaking behaviour without
# editing the script itself.
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

makedocs(
    modules = [AstroNbodySim],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaastrosim.github.io/AstroNbodySim.jl/dev/",
        assets = ["assets/alpha_small.ico"],
        analytics = "UA-153693590-1",
        highlights = ["llvm", "yaml"],
    ),
    # When true, Documenter will skip broken cross-references and
    # missing `@docs`/doctest errors instead of bailing out.
    warnonly = ("warnonly=true" in ARGS) || get(ENV, "DOCSARGS_WARN", "") == "1",
    # Documenter 1.x defaults `checkdocs` to `:all`, which makes every internal
    # docstring that is not part of the manual a hard error. Keep the pre-1.0
    # behaviour: only check what `@docs` explicitly references.
    checkdocs = :none,
    clean = false,
    sitename = "AstroNbodySim.jl",
    authors = "islent",
    linkcheck = false,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial/basic.md",
            "tutorial/mond.md",
            "tutorial/psn.md",
        ],
        "Manual" => Any[
            "manual/DirectSum.md",
            "manual/Tree.md",
            "manual/FFT.md",
            "manual/FDM.md",
            "manual/ML.md",
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
)

deploydocs(
    repo = "github.com/JuliaAstroSim/AstroNbodySim.jl.git",
    target = "build",
    devbranch = "main",
)
