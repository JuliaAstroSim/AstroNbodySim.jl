# In this example, we simulate the solar system
@info "Loading"


#=
cd("AstroNbodySim.jl/examples/07-solarsystem")
=#

# Display our solar system at present
using GLMakie
using AstroPlot
using PhysicalParticles, Unitful, UnitfulAstro
using AstroLib
using Dates
using ProgressMeter

using AstroIC
data = solarsystem(now())

unicode_scatter(data)

scene = plot_makie(data, u"AU", markersize=500)

scene, layout = plot_positionslice(solarsystem(now()), u"AU",
    xlabel = "x [AU]", yaxis = :z, ylabel = "z [AU]",
    xlims = [-22, 22], ylims = [-22, 22],
    markersize = 2.0,
    title = "Solar System at " * string(now()),
)

function display_solarsystem(startdate = now(); fps = 60.0, N = 10000, ratio = 1.0)
    last_plot_time = time()
    time_between_plot = 1.0 / fps # ms

    T = jdcnv(startdate)

    scene, layout = layoutscene()
    title = Node("Solar System at " * string(startdate))
    ax = layout[1,1] = GLMakie.Axis(
        scene,
        title = title,
        xlabel = "x [AU]",
        ylabel = "z [AU]",
        aspect = AxisAspect(1.0),
    )
    #sl1 = layout[2, 1] = Slider(scene, text = "ratio", range = 0.1:0.01:10.0, startvalue = 2.0)
    #ax2 = layout[2, 1]
    ls = labelslider!(scene, "test label: ",  0.1:0.01:10.0)
    ls.slider.value = ratio
    layout[2,1] = ls.layout

    xy = ustrip.(u"AU", pack_xy(solarsystem(T), yaxis = :z))
    GLMakie.xlims!(ax, (middle(xy[:,1]) - 22, middle(xy[:,1]) + 22))
    GLMakie.ylims!(ax, (middle(xy[:,2]) - 22, middle(xy[:,2]) + 22))

    pos = Node(xy)
    GLMakie.scatter!(ax, pos, markersize = 5.0)
    display(scene)

    @showprogress for i in 1:N
        if time() - last_plot_time > time_between_plot
            #T += to_value(sl1.value)
            T += to_value(ls.slider.value)
            xy = ustrip.(u"AU", pack_xy(solarsystem(T), yaxis = :z))
            GLMakie.xlims!(ax, (middle(xy[:,1]) - 22, middle(xy[:,1]) + 22))
            GLMakie.ylims!(ax, (middle(xy[:,2]) - 22, middle(xy[:,2]) + 22))
            @async title[] = "Solar System at " * string(daycnv(T))
            @async pos[] = xy
            last_plot_time = time()
        end
        sleep(0.1 / fps)
    end
    return scene
end

#=
display_solarsystem()
=#


# write success flag for shell
success = open("output/success", "w")
close(success)

# Clear memory
AstroNbodySim.clear()