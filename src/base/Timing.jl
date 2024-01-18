"""
$TYPEDSIGNATURES
Begin the named timer, save `time_ns()`` to it
"""
function begin_timer(sim::Simulation, name::String)
    sim.loginfo.timing[sim.loginfo.timers[name]] = 0
    return time_ns()
end

"""
$TYPEDSIGNATURES
Add time `t` to the named timer.
This is useful for procedures distributed over several functions.
"""
function add_timer(sim::Simulation, name::String, t::UInt64)
    sim.loginfo.timing[sim.loginfo.timers[name]] += t
end

"""
$TYPEDSIGNATURES
Add time `t_end - t_start` to the named timer.
This is useful for procedures distributed over several functions.
"""
function add_timer(sim::Simulation, name::String, t_start::UInt64, t_end::UInt64)
    sim.loginfo.timing[sim.loginfo.timers[name]] += t_end - t_start
end

"""
$TYPEDSIGNATURES
Save the timings to file
"""
function write_timing(sim::Simulation)
    if !(isempty(sim.loginfo.timing))
        s = ustrip(getuTime(sim.config.units), sim.timeinfo.system_time_float)
        write(sim.stream.timerio, join([s; string.(Int64.(sim.loginfo.timing))], ",") * "\n")
        flush(sim.stream.timerio)
    end
end