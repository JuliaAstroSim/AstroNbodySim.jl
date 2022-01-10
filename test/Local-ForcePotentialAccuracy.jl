### Force Accuracy
mkpathIfNotExist("output/accuracy")

TimeEnd = 0.001u"Gyr"
TimeStep = 0.00002u"Gyr"
TimeBetweenSnapshots = 0.0005u"Gyr"

pos = [PVector(u"kpc")]
append!(pos, [PVector(10.0^x, 0.0, 0.0)*1.0u"kpc" for x in -5:9])

R = ustrip.(norm.(pos))
colors = ColorSchemes.tab10.colors;

# Parameter space
Smooth = [0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0] * u"kpc"
NumParticles = [10^i for i in 1:5]
TreeOpenAngles = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]


"Relative to smoothing length"
function accuracy_force(ds, ts, pos::Vector, h::Number)
    println("h = $h")
    acc_ds = compute_force(ds, pos, h)
    acc_ts = compute_force(ts, pos, h)

    acc_delta = acc_ts - acc_ds
    ratio = (norm.(acc_delta)) ./ (norm.(acc_ds))

    return ratio
end

function accuracy_potential(ds, ts, pos::Vector, h::Number)
    println("h = $h")
    pot_ds = compute_potential(ds, pos, h)
    pot_ts = compute_potential(ts, pos, h)

    pot_delta = pot_ts - pot_ds
    ratio = (norm.(pot_delta)) ./ (norm.(pot_ds))

    return ratio
end

"Relative to particle number, at specific positions"
function accuracy_force(N::Int, pos::Vector)
    println("N = $N")
    ic = PlummerStarCluster(
        NumSamples = N,
    )
    data = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(data),
    );

    ts = Simulation(
        deepcopy(data),
        GravitySolver = Tree(),
    );

    acc_ds = compute_force(ds, pos, 0.0u"kpc")
    acc_ts = compute_force(ts, pos, 0.0u"kpc")

    acc_delta = acc_ts - acc_ds
    ratio = (norm.(acc_delta)) ./ (norm.(acc_ds))

    return ratio
end

function accuracy_potential(N::Int, pos::Vector)
    println("N = $N")
    ic = PlummerStarCluster(
        NumSamples = N,
    )
    data = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(data),
    );

    ts = Simulation(
        deepcopy(data),
        GravitySolver = Tree(),
    );

    pot_ds = compute_potential(ds, pos, 0.0u"kpc")
    pot_ts = compute_potential(ts, pos, 0.0u"kpc")

    pot_delta = pot_ts - pot_ds
    ratio = (norm.(pot_delta)) ./ (norm.(pot_ds))

    return ratio
end

"Relative to particle number, scattered"
function accuracy_force(N::Int)
    println("N = $N")
    ic = PlummerStarCluster(
        NumSamples = N,
    )
    data = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(data),
    );

    ts = Simulation(
        deepcopy(data),
        GravitySolver = Tree(),
    );

    compute_force(ds)
    compute_force(ts)

    dfds = substract_by_id(get_all_data(ds), :Acc, info = :Pos)
    dfts = substract_by_id(get_all_data(ts), :Acc, info = :Pos)

    df = diff_by_id(dfds, dfts, :Acc, info = :Pos)

    radii = ustrip.(u"kpc", norm.(df.Pos))
    ratio = (norm.(df.Acc_1 - df.Acc_2)) ./ (norm.(df.Acc_1))

    return radii, ratio
end

function accuracy_potential(N::Int)
    println("N = $N")
    ic = PlummerStarCluster(
        NumSamples = N,
    )
    data = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(data),
    );

    ts = Simulation(
        deepcopy(data),
        GravitySolver = Tree(),
    );

    compute_potential(ds)
    compute_potential(ts)

    dfds = substract_by_id(get_all_data(ds), :Potential, info = :Pos)
    dfts = substract_by_id(get_all_data(ts), :Potential, info = :Pos)

    df = diff_by_id(dfds, dfts, :Potential, info = :Pos)

    radii = ustrip.(u"kpc", norm.(df.Pos))
    ratio = (norm.(df.Potential_1 - df.Potential_2)) ./ (norm.(df.Potential_1))

    return radii, ratio
end

"Relative to opening angle"
function accuracy_force(acc_ds::Vector, pos::Vector, angle::Number, particles::StructArray)
    println("angle = $angle (rad)")

    ts = Simulation(
        deepcopy(particles),
        GravitySolver = Tree(),
        TreeOpenAngle = angle,
    );

    acc_ts = compute_force(ts, pos, 0.0u"kpc")

    acc_delta = acc_ts - acc_ds
    ratio = (norm.(acc_delta)) ./ (norm.(acc_ds))

    return ratio
end

function accuracy_potential(pot_ds::Vector, pos::Vector, angle::Number, particles::StructArray)
    println("angle = $angle (rad)")

    ts = Simulation(
        deepcopy(particles),
        GravitySolver = Tree(),
        TreeOpenAngle = angle,
    );

    pot_ts = compute_potential(ts, pos, 0.0u"kpc")

    pot_delta = pot_ts - pot_ds
    ratio = (norm.(pot_delta)) ./ (norm.(pot_ds))

    return ratio
end

function cluster_accuracy()
    @info "Relative to smoothing length"
    ic = PlummerStarCluster(
        NumSamples = 10000,
    )
    particles = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(particles),
    );

    ts = Simulation(
        deepcopy(particles),
        GravitySolver = Tree(),
    );
    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Force Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(Δf/f)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_force(ds, ts, pos, Smooth[i])), color = colors[i]) for i in eachindex(Smooth)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(Smooth), "smoothing length",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/ForceAccuracySmoothingLength.png", scene)

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Potential Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(|Δφ|/|φ|)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_potential(ds, ts, pos, Smooth[i])), color = colors[i]) for i in eachindex(Smooth)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(Smooth), "smoothing length",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/PotentialAccuracySmoothingLength.png", scene)

    @info "Relative to opening angle"
    ic = PlummerStarCluster(
        NumSamples = 10000,
    )
    particles = generate(ic, MaxRadius = 0.0u"kpc")

    ds = Simulation(
        deepcopy(particles),
    );

    acc_ds = compute_force(ds, pos, 0.0u"kpc")
    pot_ds = compute_potential(ds, pos, 0.0u"kpc")

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Force Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(Δf/f)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_force(acc_ds, pos, TreeOpenAngles[i], particles)), color = colors[i]) for i in eachindex(TreeOpenAngles)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(TreeOpenAngles), "Opening Angle",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/ForceAccuracyOpeningAngle.png", scene)

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Potential Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(|Δφ|/|φ|)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_potential(pot_ds, pos, TreeOpenAngles[i], particles)), color = colors[i]) for i in eachindex(TreeOpenAngles)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(TreeOpenAngles), "Opening Angle",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/PotentialAccuracyOpeningAngle.png", scene)
end

cluster_accuracy()

function scatter_accuracy(N::Int)
    @info "Relative to particle number"
    radii, ratio = accuracy_force(N)

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Force Error ($N particles)",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(Δf/f)",
    )
    s = Makie.scatter!(axis, log10.(radii), log10.(ratio), markersize = 5.0)
    Makie.save("./accuracy/ForceAccuracyScatter.png", scene)


    radii, ratio = accuracy_potential(N)

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Potential Error ($N particles)",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(|Δφ|/|φ|)",
    )
    s = Makie.scatter!(axis, log10.(radii), log10.(ratio), markersize = 5.0)
    Makie.save("./accuracy/PotentialAccuracyScatter.png", scene)
end

scatter_accuracy(1000)


### Single particle

# Scatter plot
function single_accuracy()
    # Add massless particles for tree construction
    ic = PlummerStarCluster(
        NumSamples = 1000,
    )
    particles = generate(ic, MaxRadius = 0.0u"kpc")        
    particles.Mass .*= 0.0

    push!(particles, Star(uAstro))
    particles.Mass[end] = 1.0e10u"Msun"

    pos = []
    for x in -1.0u"kpc":0.1u"kpc":1.0u"kpc"
        for y in -1.0u"kpc":0.1u"kpc":1.0u"kpc"
            for z in -1.0u"kpc":0.1u"kpc":1.0u"kpc"
                push!(pos, PVector(x,y,z))
            end
        end
    end
    pos = identity.(pos)
    R = ustrip.(norm.(pos))

    @info "Relative to smoothing length"

    ds = Simulation(
        deepcopy(particles),
    );

    ts = Simulation(
        deepcopy(particles),
        GravitySolver = Tree(),
    );

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Force Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(Δf/f)",
    )
    s = [Makie.scatter!(axis, log10.(R), log10.(accuracy_force(ds, ts, pos, Smooth[i])), color = colors[i]) for i in eachindex(Smooth)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(Smooth), "smoothing length",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/ForceAccuracyScatterSmoothingLength.png", scene)
    # No error, zero plot

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Potential Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(|Δφ|/|φ|)",
    )
    s = [Makie.scatter!(axis, log10.(R), log10.(accuracy_potential(ds, ts, pos, Smooth[i])), color = colors[i]) for i in eachindex(Smooth)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(Smooth), "smoothing length",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/PotentialAccuracyScatterSmoothingLength.png", scene)
    # No error, zero plot


    @info "Relative to opening angle"
    acc_ds = compute_force(ds, pos, 0.0u"kpc")
    pot_ds = compute_potential(ds, pos, 0.0u"kpc")

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Force Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(Δf/f)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_force(acc_ds, pos, TreeOpenAngles[i], particles)), color = colors[i]) for i in eachindex(TreeOpenAngles)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(TreeOpenAngles), "Opening Angle",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/ForceAccuracyScatterOpeningAngle.png", scene)
    # No error, zero plot

    scene, layout = Makie.layoutscene(resolution = (800, 450))
    axis = layout[1,1] = Makie.Axis(
        scene,
        title = "Relative Potential Error",
        xlabel = "log10(R [kpc])",
        ylabel = "log10(|Δφ|/|φ|)",
    )
    s = [Makie.lines!(axis, log10.(R), log10.(accuracy_potential(pot_ds, pos, TreeOpenAngles[i], particles)), color = colors[i]) for i in eachindex(TreeOpenAngles)]
    legend = layout[1,1] = Makie.Legend(
        scene, s, string.(TreeOpenAngles), "Opening Angle",
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
    )
    Makie.save("./accuracy/PotentialAccuracyScatterOpeningAngle.png", scene)
    # No error, zero plot
end

single_accuracy()

AstroNbodySim.clear()