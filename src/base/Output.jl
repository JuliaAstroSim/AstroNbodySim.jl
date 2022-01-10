function mkpathIfNotExist(dir)
    if !isdir(dir)
        mkpath(dir)
    end
end

function need_to_interrupt(OutputDir::String)
    if isfile(joinpath(OutputDir, "stop"))
        return true
    else
        return false
    end
end

function count_gadget_types(all_data::StructArray) where T<:AbstractParticle
    Counts = MVector{6,Int32}([0,0,0,0,0,0])
    for p in all_data.Collection
        Counts[Int(p)] += 1
    end
    return Counts
end

function HeaderGadget2(sim::Simulation, all_data::StructArray) where T<:AbstractParticle
    Counts = count_gadget_types(all_data)
    return HeaderGadget2(
        Counts,
        MVector{6,Float64}([0.0,0.0,0.0,0.0,0.0,0.0]),
        ustrip(getuTime(getunits(sim)), sim.timeinfo.system_time_float),
        sim.timeinfo.redshift, 0, 0,
        Counts,
        0,
        1,

        0.0, 0.3, 0.7, 0.71, 0, 0,
        MVector{6,UInt32}([0,0,0,0,0,0]), 0,
        @MVector zeros(UInt8, 60)
    )
end

function output(sim::Simulation, filename::String, OutputType::jld2)
    all_data = get_all_data(sim)
    if sim.config.output.pot
        compute_potential(sim)
    end
    write_jld(filename * ".jld2", all_data)
end

function output(sim::Simulation, filename::String, OutputType::gadget2; kw...)
    all_data = get_all_data(sim)
    header = HeaderGadget2(sim, all_data)

    acc = sim.config.output.acc
    pot = sim.config.output.pot
    if pot
        compute_potential(sim)
    end

    write_gadget2(filename * ".gadget2", header, Measurements.value.(all_data), sim.config.units;
        format2 = sim.config.output.format2,
        acc, pot,
        kw...
    )
end

function output(sim::Simulation, OutputType::AbstractOutputType)
    OutputDir = sim.config.output.dir
    if !isdir(OutputDir)
        @warn "Output Directory does not exist! Making a new one" mkpath(OutputDir)
    end
    filename = string(OutputDir, "/snapshot_", @sprintf "%04d" sim.outputinfo.snapshotcount)

    output(sim, filename, OutputType)
    
    sim.outputinfo.snapshotcount += 1
end

function outputparallel(sim::Simulation, OutputType::AbstractOutputType)
end