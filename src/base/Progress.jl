function setup_progress(sim::Simulation, TimestepConfig::AdaptiveTimestep, color = :light_blue; kw...)
    sim.visinfo.progress = Progress(TimestepConfig.TimeBase; #=showspeed=true,=# color)
end

function setup_progress(sim::Simulation, TimestepConfig::ConstantTimestep, color = :light_blue; kw...)
    totalstep = ceil(Int64, (sim.config.time.End - sim.config.time.Begin) / TimestepConfig.dt)
    sim.visinfo.progress = Progress(totalstep; #=showspeed=true,=# color)
end

function update_progress(sim::Simulation, TimestepConfig::AdaptiveTimestep; kw...)
    ProgressMeter.update!(sim.visinfo.progress, sim.timeinfo.system_time_int; kw...)
end

function update_progress(sim::Simulation, TimestepConfig::ConstantTimestep; kw...)
    time = floor(Int64, (sim.timeinfo.system_time_float - sim.config.time.Begin) / TimestepConfig.dt)
    ProgressMeter.update!(sim.visinfo.progress, time; kw...)
end