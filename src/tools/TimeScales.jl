"""
    typicalvelocity(G, M, R)

v = sqrt(G*M/R)

Typical velocity of a spherical system.
"""
function typicalvelocity(M, R, G)
    return sqrt(G*M/R)
end

"""
    meandensity(M, R)

ρ = 3*M/(4πR^3)

Mean density of a spherical system.
"""
function meandensity(M, R)
    
end

"""
    crosstime(R, v)

t_cross = R / v

The typical time needed to cross the system.
"""
function crosstime(R, v)
    return R/v
end

"""
    hubbletime(H0)

t_hubble = 1/H0

The age of the universe.
"""
function hubbletime(H0)
    return 1/H0
end

"""
    relaxtime(R, v, N)

t_relax = N/(10 lnN) t_cross

Relaxation time: The time over which the change in kinetic energy due to the long-range collisions
has accumulated to a value that is comparable to the intrinsic kinetic energy of the particle.
"""
function relaxtime(R, v, N)
    return N / (10 * log10(N)) * crosstime(R, v)
end

"""
interactiontime(R, v, N)

The typical time between two short-range interactions that cause a change in kinetic energy
comparable to the intrinsic kinetic energy of the particle.

t_interaction = N * t_cross
"""
function interactiontime(R, v, N)
    return N * crosstime(R, v)
end

"""
    dynamicaltime(ρ, G)

t_dyn = sqrt(3 * π / (16 * G * ρ))

The time required to travel halfway across the system.
"""
function dynamicaltime(ρ, G)
    return sqrt(3 * π / (16 * G * ρ))
end

"""
    freefalltime(ρ, G)

t_ff = sqrt(3 * π / (32 * G * ρ))

The time it takes a sphere with zero pressure to collapse to a point.
"""
function freefalltime(ρ, G)
    return sqrt(3 * π / (32 * G * ρ))
end

"""
    orbitaltime(ρ, G)

t_orb = sqrt(3*π/(G*ρ))

The time it takes to complete a (circular) orbit.
"""
function orbitaltime(ρ, G)
    return sqrt(3*π/(G*ρ))
end
