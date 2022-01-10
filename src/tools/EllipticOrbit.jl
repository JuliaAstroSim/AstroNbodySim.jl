"""
$(TYPEDSIGNATURES)

Length of semi-major axis of elliptic orbit

- `G`: gravitational constant
- `M`: mass of central object
- `r`: orbit radius at specific time
- `v`: velocity at radius `r`
"""
function ellipticSemiMajor(G::Number, M::Number, r::Number, v::Number)
    return inv(2 / r - v^2 / (G * M))
end

"""
$(TYPEDSIGNATURES)

Orbital period of elliptic orbit

- `G`: gravitational constant
- `M`: mass of central object
- `a`: length of semi-major axis
"""
function ellipticPeriod(G::Number, M::Number, a::Number)
    return 2 * pi * sqrt(a^3 / (G * M))
end

"""
$(TYPEDSIGNATURES)

Eccentricity of elliptic orbit

- `G`: gravitational constant
- `M`: mass of central object
- `v`: velocity at radius `r`
- `h`: specific relative angular momentum at radius `r`
"""
function eccentricity(G::Number, M::Number, v::Number, h::Number)
    return 1.0 - (v * h) / (G * M)
end

"""
$(TYPEDSIGNATURES)

Eccentricity vector of elliptic orbit

- `G`: gravitational constant
- `M`: mass of central object
- `r`: position at specific time
- `v`: velocity vector at position `r`
"""
function eccentricity(G::Number, M::Number, r::AbstractPoint, v::AbstractPoint)
    return normalize(r) - cross(v, cross(r, v)) / (G * M)
end
