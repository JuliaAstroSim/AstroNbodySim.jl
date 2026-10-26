# FDM Poisson solver

The `FDM` (finite-difference) solver discretises the Poisson equation on a
Cartesian mesh and solves the resulting linear system directly. It is the
natural choice for one- and two-dimensional verification problems and for
boundary conditions that are not periodic. The matrix construction and the
linear solve are provided by
[`PhysicalFDM.jl`](https://github.com/JuliaAstroSim/PhysicalFDM.jl); the
high-level entry point is `fdm_poisson`.

## Discretisation

With the second-order central-difference Laplacian,

```math
\frac{\phi_{i-1} - 2\phi_i + \phi_{i+1}}{\Delta^2} = 4\pi G \rho_i,
```

the solver assembles ``A\mathbf{\phi} = \mathbf{b}`` with
``\mathbf{b} = 4\pi G\,\rho``. After the potential is obtained,
[`compute_acc`](@ref) evaluates the acceleration by central differences of
``\phi``.

The matrix has ``(N_x + 2N_G + 1)(N_y + 2N_G + 1)(N_z + 2N_G + 1)``
unknowns. Dense
storage grows as ``\mathcal{O}(N^{2d})`` in ``d`` dimensions; even sparse
storage suffers from fill-in during the LU factorisation, and the solve
remains expensive. This is why FDM is meant for 1D and 2D
verification problems, while 3D galaxy simulations should use the FFT solver.
The current implementation rebuilds and factorises the matrix
at every time step; the factorisation is not reused.

## Boundary conditions

| Boundary | Support |
| -------- | ------- |
| `Periodic()` | cyclic stencil |
| `Dirichlet()` | fixed (zero) boundary values |
| `Vacuum()` | open boundary, boundary potential set explicitly on the outer layer |

Neumann and Robin conditions are not implemented in `fdm_poisson`. The
`sparse` keyword is honoured on the CPU; `sparse = true` combined with
`device = GPU()` is rejected by `check_compatibility`.

## Configuration

```julia
pm = Simulation(particles;
                GravitySolver     = FDM(),
                device            = CPU(),
                Nx = Ny = Nz      = 32,
                xMin = -50u"kpc", xMax = +50u"kpc",
                NG                = 1,
                BoundaryCondition = Dirichlet(),
                sparse            = true)
run(pm)
```

For 2D and 3D runs the same matrix machinery is used; the mesh dimension is
inferred from the simulation (`dim = 3` for particle data, or `1`/`2` for
explicitly constructed `MeshCartesianStatic` test meshes).

## API

The solver and the finite-difference operators live in `PhysicalFDM.jl`:

- `fdm_poisson(mesh, Val(dim), G, mode, device, boundary, sparse)` — assemble
  and solve the Poisson system.
- `compute_acc(mesh)` — acceleration by central differences of `mesh.phi`
  (implemented in `AstroNbodySim.jl`).
- `diff_mat` — finite-difference matrix for a given stencil and boundary.
- `laplace_conv_op` — convolution form of the Laplacian used by the
  finite-difference operators.
