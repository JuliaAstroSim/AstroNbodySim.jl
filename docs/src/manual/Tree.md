# Tree gravity solver

The `Tree` solver is a Barnes–Hut octree
([Barnes & Hut 1986, Nature 324, 446](https://doi.org/10.1038/324446a0))
that reduces the pairwise force calculation to ``\mathcal{O}(N \log N)``.
The octree itself is implemented in
[`PhysicalTrees.jl`](https://github.com/JuliaAstroSim/PhysicalTrees.jl);
the force walk and the distributed-memory bookkeeping follow the Gadget-2
framework (Springel, Yoshida & White 2001; Springel 2005).

## When to use it

- ``N \gtrsim 10^5`` particles, where the ``N^2`` DirectSum cost becomes
  prohibitive.
- Runs that are already distributed over several processes: the Peano–Hilbert
  top-tree and pseudo-node mechanism are designed for that case.

## Tree construction

1. Particles are sorted by their Peano–Hilbert key inside a domain that is
   slightly enlarged by `ExtentMargin` (default `1.001`).
2. The sorted sequence is cut into top-level nodes (`TopnodeFactor`,
   default `20`). Each top-level node is a contiguous block of the Peano
   curve and becomes the unit of domain decomposition.
3. Inside each domain the particles are inserted recursively into a cubic
   octree. Every node stores its mass, centre of mass and side length; the
   force walk uses these as a monopole approximation.
4. Nodes owned by another process are represented by **pseudo-nodes**
   (indices above `MaxTreenode`). When the walk of a local target reaches a
   pseudo-node, the target is exported to the owning process; that process
   then applies the opening criterion to its own tree. Remote tree data are
   therefore never fetched by the local process.

The tree is rebuilt from scratch after every drift (`rebuild` in
`src/tree/timestep.jl`); incremental tree updates are not implemented.

## Opening criteria

Two criteria are available. The default is the geometric Barnes–Hut
criterion; the Gadget-2 relative criterion is selected by setting
`TreeOpenAngle <= 0`.

**Geometric criterion (`TreeOpenAngle > 0`, default `0.1` rad).**
A node of side length ``l`` at distance ``r`` from the target is accepted as a
single monopole contribution when

```math
\frac{l}{r} < \theta, \qquad \theta = \texttt{TreeOpenAngle},
```

and is opened otherwise. The implementation tests the equivalent form
``l^2 > r^2 \theta^2``.

**Relative criterion (`TreeOpenAngle <= 0`).**
A node is opened when

```math
\frac{G M}{r^2}\left(\frac{l}{r}\right)^2
> \alpha\,\lVert \mathbf{a}_{\mathrm{old}} \rVert,
\qquad \alpha = \texttt{ErrTolAcc},
```

where ``M`` is the node mass, ``l`` its side length, and
``\mathbf{a}_{\mathrm{old}}`` the acceleration of the target from the previous
step. In addition, a node is always opened when the target lies inside
``0.6\,l`` of the node centre. `ErrTolAcc` defaults to `0.025`.

Some formulations write this criterion with a single factor ``l/r``; the
implementation and the Gadget-2 reference use the squared form
``(l/r)^2`` shown above.

For either criterion, a node that is accepted contributes

```math
\mathbf{a}_i \mathrel{+}= G M_{\mathrm{node}}\,
\frac{\mathbf{r}_{\mathrm{node}} - \mathbf{r}_i}{h_i^3}\,
g\!\left(\frac{\lVert \mathbf{r}_{\mathrm{node}} - \mathbf{r}_i \rVert}{h_i}\right),
```

with the same Gadget-2 cubic-spline kernel ``g`` and per-collection softening
length ``h_i`` used by the [DirectSum solver](@ref "DirectSum gravity solver").
As ``\theta \to 0`` (or ``\alpha \to 0``) the tree walk degenerates to
DirectSum.

## Configuration

The most useful keywords accepted by `Simulation` are listed below. The
remaining allocation parameters (`MaxTopnode`, `TopnodeFactor`,
`ToptreeAllocFactor`, `TreeAllocFactor`, `ExtentMargin`, `PeanoBits3D`,
`PeanoBits2D`, `epsilon`) are forwarded to `PhysicalTrees.OctreeConfig`.

| Keyword | Default | Description |
| ------- | ------- | ----------- |
| `TreeOpenAngle` | `0.1` rad | Barnes–Hut opening angle. Set to `0` or negative to switch to the relative criterion. |
| `ErrTolAcc` | `0.025` | Accuracy parameter ``\alpha`` of the relative criterion. |
| `MaxTreenode` | `10N` (`N > 1000`), else `10000` | Capacity of the local octree. |
| `ExtentMargin` | `1.001` | Fractional padding of the tree domain. |
| `PeanoBits3D` | `21` | Peano–Hilbert key length per axis (3D). |

```julia
ts = Simulation(particles;
                GravitySolver       = Tree(),
                TreeOpenAngle       = 0.1,   # rad
                ErrTolAcc           = 0.025,
                ForceSofteningTable = fill(0.1u"kpc", 6))
run(ts)
```

## Distributed memory

The top-tree is the communication skeleton. Each process walks its local
octree for its own targets; whenever the walk reaches a pseudo-node, the
target position is appended to `sim.buffer.gravToEval[pid]` for the owning
process. The remote process evaluates those targets against its local tree and
returns `GravResult` entries, which `combine_force` adds to `data.Acc`.
Communication uses `Distributed.jl` (TCP/IP); an MPI backend is not
implemented.

Time integration and timestep control are shared with the DirectSum solver:
Gadget-2 leapfrog (kick–drift–kick), with `TimeStep == 0` selecting the
adaptive step. See [DirectSum gravity solver](@ref "DirectSum gravity solver")
for the equations.

## References

- Barnes, J. & Hut, P. 1986, Nature, 324, 446.
- Springel, V., Yoshida, N. & White, S. D. M. 2001, New Astronomy, 6, 79.
- Springel, V. 2005, Nature, 435, 629.

## API

```@docs
TreeSimConfig
compute_force
softlen
```

See [`PhysicalTrees.jl`](https://juliaastrosim.github.io/PhysicalTrees.jl/dev)
for the low-level octree, Peano-key and domain-decomposition API.
