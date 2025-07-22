"""
    apply_background_force(f::Function, data::StructArray)

Apply `f` to all particles in data by calling `data.Acc[i] += f(data[i], system_time_float)`.
`f` must return a `PVector` in the unit of physical acceleration.
"""
Base.@propagate_inbounds function apply_background_force(f::Function, data::StructArray, system_time_int::Int, system_time_float)
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds te = data.Ti_endstep[i]
            if te == system_time_int
                @inbounds data.Acc[i] += f(data[i], system_time_float)
            end
        end
    end
end

"""
    compute_OldAcc(data::StructArray)

`data.OldAcc[i] = norm(data.Acc[i])`
"""
Base.@propagate_inbounds function compute_OldAcc(data::StructArray)
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            data.OldAcc[i] = norm(data.Acc[i])
        end
    end
end

function postprocessing_force(sim::Simulation, ::Gravity, ::CPU)
    #! Do not change the order!
    data = get_local_data(sim)

    for f in sim.bgforce
        apply_background_force(f, data, sim.timeinfo.system_time_int, sim.timeinfo.system_time_float)
    end

    # MOND - Milgrom 1983
    if sim.config.grav.model isa MOND1983Milgrom
        mond_Milgrom1983(sim, data, sim.timeinfo.system_time_int)
    end

    compute_OldAcc(data)
end