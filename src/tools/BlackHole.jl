"""
    r_g(G, M, c)

r_g = 2 * G * M / c^2

where `G` is the gravity constant, `M` is the mass of the central object, `c` is light speed.

Schwarzchild radius.
"""
function r_g(G, M, c)
    return 2 * G * M / c^2
end

radius_gravity = r_g
radius_schwarzschild = r_g


"""
    pseudoNewtonianPotential(G, M, c, R)

pot = - G * M / (R - r_g(G, M, c))

Pseudo-Newtonian potential around a compat central object. Diverge at gravity radius `r_g`.
"""
function pseudoNewtonianPotential(G, M, c, R)
    return - G * M / (R - r_g(G, M, c))
end

"""
pseudoNewtonianPotential(G, M, rg)

pot = - G * M / (R - rg)

Pseudo-Newtonian potential around a compat central object. Diverge at gravity radius `r_g`.
"""
function pseudoNewtonianPotential(G, M, rg)
    return - G * M / (R - rg)
end


"""
    pseudoNewtonianAcc(G, M, c, R, n::AbstractPoint)

acc = - G * M / (R - r_g(G, M, c))^2 * n

Pseudo-Newtonian acceleration around a compat central object. Diverge at gravity radius `r_g`.
"""
function pseudoNewtonianAcc(G, M, c, R, n::AbstractPoint)
    return - G * M / (R - r_g(G, M, c))^2 * n
end

"""
    pseudoNewtonianAcc(G, M, c, R, n::AbstractPoint)

acc = - G * M / (R - rg)^2 * n

Pseudo-Newtonian acceleration around a compat central object. Diverge at gravity radius `r_g`.
"""
function pseudoNewtonianAcc(G, M, rg, n::AbstractPoint)
    return - G * M / (R - rg)^2 * n
end