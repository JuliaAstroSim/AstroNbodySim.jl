csvstring(d::Any) = string(d)

function csvstring(d::AbstractPoint)
    return string("\"", d, "\"")
end

function write_analysis(sim::Simulation)
    if !isempty(sim.loginfo.analysers)
        s = string(ustrip(getuTime(sim.config.units), sim.timeinfo.system_time_float), ",")
        analysis = []
        push!(analysis, sim.timeinfo.system_time_float)
        for key in keys(sim.loginfo.analysers)
            func = sim.loginfo.analysers[key]
            d = func(sim)
            ds = csvstring(ustrip(d))
            s = string(s, ds, ",")
            push!(analysis, d)
        end
        write(sim.stream.analyserio, chop(s) * "\n")
        flush(sim.stream.analyserio)
        return analysis
    end
end