function plot_training(filename::AbstractString;
    xlabel = "iter",
    ylabel = "loss",
    resolution = (1600, 900),
    title = "Training loss curve",
    colors = ColorSchemes.tab10.colors,
    yscale = log10,
    kw...
)
    df = DataFrame(CSV.File(filename))
    fig = Figure(; resolution)
    ax = Makie.Axis(fig[1,1]; xlabel, ylabel, title, yscale)
    p1 = Makie.lines!(ax, df.iter, df.MSEloss; kw...)
    p2 = Makie.lines!(ax, df.iter, df.Maxloss; kw...)
    Legend(fig[1,2], [p1, p2], ["MSEloss", "Maxloss"])
    Makie.save(filename*".png", fig)

    return fig
end