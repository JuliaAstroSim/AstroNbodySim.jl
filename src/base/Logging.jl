function setuploggers(sim::Simulation)
    stream = sim.stream
    loginfo = sim.loginfo
    config = sim.config

    # stream.loggingio = open(joinpath(config.output.dir, "sim.log"), "w")

    if !isempty(loginfo.timing)
        stream.timerio = open(joinpath(config.output.dir, "timing.csv"), "w")
        write(stream.timerio, join(["time"; string.(keys(loginfo.timers))...], ",") * "\n")
    end

    if !isempty(loginfo.analysers)
        stream.analyserio = open(joinpath(config.output.dir, "analysis.csv"), "w")
        write(stream.analyserio, join(["time"; string.(keys(loginfo.analysers))], ",") * "\n")
    end

    if sim.config.loggingmode == ProgressMode()
        setup_progress(sim, sim.config.time.step)
    end
end

function restartloggers(sim::Simulation)
    stream = sim.stream
    loginfo = sim.loginfo
    config = sim.config

    # sim.loggingio = open(joinpath(config.output.dir, "sim.log"), "a")
    
    if !isempty(loginfo.timing)
        stream.timerio = open(joinpath(config.output.dir, "timing.csv"), "a")
    end

    if !isempty(loginfo.analysers)
        stream.analyserio = open(joinpath(config.output.dir, "analysis.csv"), "a")
    end

    if sim.config.loggingmode == ProgressMode()
        setup_progress(sim, sim.TimestepConfig)
    end

    #TODO restart progress
end

function closeloggers(sim::Simulation)
    stream = sim.stream
    loginfo = sim.loginfo

    #! Before saving, assign io hundlers (which are pointers) to nothing.
    # close(stream.loggingio)
    close(stream.timerio)
    close(stream.analyserio)

    # Progress is also a pointer type
    sim.visinfo.progress = Progress(0)
end

# function write_log(sim::Simulation, args...; color::Union{Symbol,Int}=:normal)
#     write(sim.stream.loggingio, args...)
#     flush(sim.stream.loggingio)

#     if sim.config.loggingmode == NormalMode()
#         printstyled(args...; color = color)
#     end
# end

#TODO
function error_log(sim::Simulation, args...)
    e = open(joinpath(sim.config.output.dir, string("error.step_", sim.timeinfo.stepcount, ".log")))
    write(e, args...)
    close(e)

    error(args...)
end