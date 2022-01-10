"""
    nu(y::Number, n::Number = 2)

Generalized MOND (MOdified Newtonian Dynamics) interpolation function. See also `nu1`, `nu2`
"""
function nu(y::Number, n::Number = 2)
    return (0.5 * (1.0 + sqrt(1.0 + 4.0 * float(y)^-n))) ^ (1.0 / n)
end

"""
    nu1(y::Number)

MOND (MOdified Newtonian Dynamics) interpolation function with index = 1. See also `nu`, `nu2`
"""
function nu1(y::Number)
    return 0.5 * (1.0 + sqrt(1.0 + 4.0 / y))
end

"""
    nu1(y::Number)

MOND (MOdified Newtonian Dynamics) interpolation function with index = 2. See also `nu`, `nu1`
"""
function nu2(y::Number)
    return sqrt(0.5 * (1.0 + sqrt(1.0 + 4.0 / y^2)))
end

"""
    mond_Milgrom1983(sim::Simulation, data::StructArray)

Apply Milgrom 1983 formula of MOND (MOdified Newtonian Dynamics) to accelerations
"""
function mond_Milgrom1983(sim::Simulation, data::StructArray)
    nuIndex = sim.config.grav.MOND_nuIndex
    ACC0 = sim.config.constants.ACC0

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds NewtonAcc = data[i].Acc
            @inbounds data.Acc[i] = NewtonAcc * nu(norm(NewtonAcc) / ACC0, nuIndex)
        end
    end
end