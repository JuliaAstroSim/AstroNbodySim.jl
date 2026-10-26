# FFT Poisson solver

The `FFT` solver treats gravity as a Poisson problem on a Cartesian mesh and
computes the convolution with the Green's function by fast Fourier transform.
The implementation is provided by
[`PhysicalFFT.jl`](https://github.com/JuliaAstroSim/PhysicalFFT.jl) and is
exposed here through `fft_poisson` / `fft_poisson!`.

## Discretisation

In the continuum the potential satisfies

```math
\nabla^2 \phi = 4\pi G \rho,
\qquad
\widehat{\phi}(\mathbf{k}) =
-\frac{4\pi G\,\widehat{\rho}(\mathbf{k})}{\lVert \mathbf{k} \rVert^2}.
```

`PhysicalFFT.jl` uses the *discrete* Fourier symbol of the second-order
Laplacian, so on a uniform grid with spacing ``\Delta_x, \Delta_y, \Delta_z``
the denominator of the periodic solver is

```math
\frac{2}{\Delta_x^2}\left(\cos k_x\Delta_x - 1\right)
+\frac{2}{\Delta_y^2}\left(\cos k_y\Delta_y - 1\right)
+\frac{2}{\Delta_z^2}\left(\cos k_z\Delta_z - 1\right).
```

The ``\mathbf{k} = 0`` mode is explicitly set to zero, which is equivalent to
solving for ``\rho - \langle \rho \rangle``; this is the same convention as in
the wave-DM Poisson step of the WaveDM.jl framework.

## Boundary conditions

| Boundary | Implementation | Dimensions |
| -------- | -------------- | ---------- |
| `Periodic()` | FFT with the discrete Laplacian symbol above | 1D / 2D / 3D |
| `Dirichlet()` | homogeneous Dirichlet, discrete sine transform (`RODFT00`) | 1D / 2D / 3D |
| `Vacuum()` | James (1977) screening-charge correction | 3D, CPU |

Only the periodic path has GPU implementations; the Dirichlet and vacuum
solvers are CPU-only in the current `PhysicalFFT.jl`.

The vacuum solver first solves the homogeneous Dirichlet problem, computes
the "screening charges" on the six faces of the box, and convolves them with
the Green's function by direct summation. Its cost is ``\mathcal{O}(N^4)``
for an ``N^3`` mesh, so it is intended for verification and low-resolution
runs. For production isolated-halo runs the periodic
solver is recommended, with a box at least eight times the halo scale radius,
so that negligible mass crosses the boundary.

`Simulation` defaults to `BoundaryCondition = Vacuum()`; pass
`BoundaryCondition = Periodic()` explicitly for the FFT production setup.
Note that the vacuum FFT path is only implemented in three dimensions.

## Configuration

```julia
pm = Simulation(particles;
                GravitySolver     = FFT(),
                device            = CPU(),
                Nx = Ny = Nz      = 64,
                xMin = -50u"kpc", xMax = +50u"kpc",
                NG                = 1,
                BoundaryCondition = Periodic())
run(pm)
```

`NG` is the number of ghost cells per side. The mesh extent is enlarged by
`EnlargeMesh` (default `2.01`) when it is derived from the particle
distribution. `device = GPU()` selects the CUDA FFT path for the periodic
solver.

## Parallel execution

FFT operations parallelise over shared-memory threads and are the recommended
target for GPU execution:

```julia
using FFTW
FFTW.set_num_threads(Threads.nthreads())
FFTW.set_provider!("mkl")   # optional, if MKL is available
```

Distributed-memory FFT (MPI) is not implemented in this release.

See `examples/06-GalaxyCollision` for a disk–disk collision. Mesh quantities
can be sampled back to particle positions with `PhysicalMeshes.mesh2particle`.

## References

- James, R. A. 1977, Journal of Computational Physics, 25, 71 (vacuum boundary correction).

## API

The implementation lives in `PhysicalFFT.jl`:

- `fft_poisson(Δ, Len, rho, boundary, device)` — solve the Poisson equation
  for a raw density array with spacing `Δ` and array size `Len`.
- `fft_poisson!(mesh, rho, boundary)` — in-place version writing to
  `mesh.phi`.
- `fft_poisson(mesh, G)` — solve for `mesh.rho` with gravitational constant
  `G`, using `mesh.config.boundary`.
