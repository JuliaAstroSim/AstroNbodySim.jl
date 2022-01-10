csvstring(d::Any) = string(d)

function csvstring(d::AbstractPoint)
    return string("\"", d, "\"")
end

function write_analysis(sim::Simulation)
    if !isempty(sim.loginfo.analysers)
        s = string(ustrip(getuTime(sim.config.units), sim.timeinfo.system_time_float), ",")
        for key in keys(sim.loginfo.analysers)
            func = sim.loginfo.analysers[key]
            d = ustrip(func(sim))
            ds = csvstring(d)
            s = string(s, ds, ",")
        end
        write(sim.stream.analyserio, chop(s) * "\n")
        flush(sim.stream.analyserio)
    end
end