# Tutorial: pseudo-Newtonian potentials

`AstroNbodySim.jl` provides a small set of helpers for integrating orbits near
a compact object with a pseudo-Newtonian (PSN) potential. They are ordinary
functions, not a `Simulation` keyword: you call them directly or wrap them in
a background-force function (see below).

## Definitions

The gravitational radius used here is the Schwarzschild radius

```math
r_g = \frac{2GM}{c^2},
```

and the potential and acceleration are

```math
\Phi_{\mathrm{PSN}}(R) = -\frac{GM}{R-r_g},
\qquad
\mathbf{a}_{\mathrm{PSN}}(R) = -\frac{GM}{(R-r_g)^2}\,\hat{\mathbf{n}}.
```

This is the Paczyński–Wiita form, which reproduces the Schwarzschild
periastron precession to leading order. The potential diverges at
``R = r_g``; only orbits with ``R > r_g`` are meaningful.

## Available functions

| Function | Signature |
| -------- | --------- |
| [`r_g`](@ref) | `r_g(G, M, c)` |
| `radius_gravity`, `radius_schwarzschild` | aliases of `r_g` |
| [`pseudoNewtonianPotential`](@ref) | `pseudoNewtonianPotential(G, M, c, R)` |
| [`pseudoNewtonianAcc`](@ref) | `pseudoNewtonianAcc(G, M, c, R, n)` |

Here `n` is a unit vector (for example `normalize(p.Pos)`) pointing away from
the compact object; the returned acceleration points inward.

## Using it as a background force

```julia
using AstroNbodySim, PhysicalParticles
using Unitful, UnitfulAstro
astro()

const G = AstroNbodySim.Constant(uAstro).G
const c = AstroNbodySim.Constant(uAstro).c
const M = 1.0e8u"Msun"

function psn_force(p::AbstractParticle, t)
    R = norm(p.Pos)
    n = normalize(p.Pos)
    return pseudoNewtonianAcc(G, M, c, R, n)
end

sim = Simulation(particles;
                 bgforce  = Function[psn_force],
                 TimeEnd  = 0.1u"Gyr",
                 TimeStep = 0.0u"Gyr")
run(sim)
```

The background force is added to `data.Acc` at every force evaluation; the
required signature is `f(particle, time) -> acceleration`.

## API

`r_g`, `pseudoNewtonianPotential` and `pseudoNewtonianAcc` are documented on
the Library → Methods page.
