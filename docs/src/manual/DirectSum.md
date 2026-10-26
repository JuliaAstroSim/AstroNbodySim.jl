# DirectSum gravity solver

`DirectSum` evaluates the gravitational acceleration of every particle by
summing over all source particles — an ``\mathcal{O}(N^2)`` operation that is
free of multipole truncation error and straightforwardly parallel.

## When to use it

- **Small and medium problems.** Published benchmarks use DirectSum for
  ``N \lesssim 10^5`` particles; on a single CPU core the practical limit is
  closer to ``N \sim 10^4`` per step, while multi-threading, distributed-memory
  execution or a GPU extend the useful range.
- **Reference solutions.** DirectSum is the ground truth against which the
  Tree, FFT and FDM solvers are compared in the tests and in
  `examples/05-plummer`.
- **Kernel prototyping.** New force kernels (tidal tensors, post-Newtonian
  corrections, custom background fields) are easiest to validate here first.

## Algorithm

For a target particle ``i`` with softening length ``h_i`` and sources ``j``,

```math
\mathbf{a}_i = G \sum_{j \neq i} m_j \,
               \frac{\mathbf{r}_j - \mathbf{r}_i}{h_i^3} \,
               g\!\left(\frac{\lVert \mathbf{r}_j - \mathbf{r}_i \rVert}{h_i}\right),
```

where ``g(u)`` is the cubic-spline softening kernel used by Gadget-2
(Springel, Yoshida & White 2001; Springel 2005):

```math
g(u) =
\begin{cases}
\dfrac{32}{3} - 38.4\,u^2 + 32\,u^3, & 0 \le u < \tfrac{1}{2}, \\
\dfrac{64}{3} - 48\,u + 38.4\,u^2 - \dfrac{32}{3}\,u^3
 - \dfrac{1}{15\,u^3}, & \tfrac{1}{2} \le u < 1, \\
\dfrac{1}{u^3}, & u \ge 1.
\end{cases}
```

The kernel is ``\mathcal{C}^2`` at ``u = 1/2`` and ``u = 1``. For ``u \ge 1``
the expression reduces to the usual Newtonian acceleration
``G m_j (\mathbf{r}_j - \mathbf{r}_i)/r^3``.

The implementation is in
[`src/directsum/gravity.jl`](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/src/directsum/gravity.jl).
Note that the helper `compute_force_pairwise` uses a *Plummer*
softening ``(\lVert\Delta\mathbf{r}\rVert^2 + h^2)^{-3/2}``; it is not the
kernel used by the main solver loop above. Some references describe the force
softening generically as Plummer softening, while the production kernel in
`compute_force_from_particle` is the Gadget-2 cubic spline given above.

## Softening length

The softening length is selected per particle `Collection` from the
six-entry `ForceSofteningTable` (gas, halo, disk, bulge, star, black hole).
`SimConfig` defaults to `1.0e-4` length units (kpc with `uAstro`) for every
collection; `GravityConfig` uses `0.01u"kpc"`.

A useful starting estimate is the mean inter-particle spacing
(the choice adopted by `suggest_softlen`),

```math
h \sim \left(\frac{V}{N}\right)^{1/3},
```

which is exactly [`suggest_softlen`](@ref). Set it explicitly with
[`set_softlen!`](@ref) or [`suggest_softlen!`](@ref). The value must be
strictly positive and is used both in the force kernel and in the adaptive
timestep criterion.

## Multi-threading, distributed memory and GPU

The three parallel modes are implemented as follows.

| Mode | Implementation |
| ---- | -------------- |
| Multi-threading | The local target list is split into `Threads.nthreads()` blocks and evaluated with `Threads.@threads`. |
| Distributed memory | Each process computes the contribution of its local sources to targets owned by other processes. Targets are staged in `sim.buffer.gravToEval`, results return in `sim.buffer.gravResult`, and `combine_force` adds them to `data.Acc`. Communication uses `Distributed.jl` (TCP/IP), not MPI. |
| GPU | A single all-pairs CUDA kernel assigns one thread per target particle and loops over all sources (`GPU_Algorithm = AllPairs()`). `Tiled` is declared as a trait but has no implementation yet. |

Because every pair contribution is independent, the summation order only
changes round-off at the last bit, which is what makes the distributed
decomposition exact up to floating-point associativity.

## Time integration

Particle solvers share the Gadget-2-style leapfrog (kick–drift–kick) scheme:
a half kick, a drift, a force evaluation and a second half kick. The default
is `TimeIntegrationAlgorithm = Leapfrog()`; `Euler()` is also available.

- `TimeStep > 0` selects a constant step.
- `TimeStep == 0` (the `SimConfig` default) selects the adaptive step. The
  step is quantised on an integer timeline with `TimeBase = 1 << 30`, and the
  raw step is

  ```math
  \Delta t = \sqrt{\frac{2\,\eta\,h_i}{\lVert \mathbf{a}_i \rVert}},
  ```

  with `ErrTolTimestep = 0.025` and `MinStep`/`MaxStep` as bounds.

## Static background force

A common use of DirectSum and Tree is to freeze a particle distribution and
evaluate its force at arbitrary positions, for example to build the static
baryon background used in `WaveDM.jl`'s `:particles_static` mode:

```julia
baryon_sim = Simulation(baryon_particles; GravitySolver = Tree())
acc = compute_force(baryon_sim, positions, 0.1u"kpc")  # Vector of PVector
```

The position-based form does not advance the particles; it only evaluates the
softened force at the supplied points.

## Typical configuration

```julia
ds = Simulation(particles;
                GravitySolver = DirectSum(),
                device        = CPU(),
                ForceSofteningTable = fill(0.1u"kpc", 6),
                TimeEnd       = 0.1u"Gyr",
                TimeStep      = 0.0u"Gyr")  # 0 → adaptive
run(ds)
```

See `examples/05-plummer.md` for Plummer-cluster runs and the
[regression tests](https://github.com/JuliaAstroSim/AstroNbodySim.jl/tree/main/test)
for quantitative cross-solver comparisons.

## References

- Springel, V., Yoshida, N. & White, S. D. M. 2001, New Astronomy, 6, 79.
- Springel, V. 2005, Nature, 435, 629.
