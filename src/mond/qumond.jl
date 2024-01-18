"""
    QUMOND_PDM_density(m::MeshCartesianStatic, ACC0::Number)

Compute ρ_PDM on the RHS (right hand side) of QUMOND (QUasi-linear MOdified Newtonian Dynamics).
Return ρ_PDM
"""
function QUMOND_PDM_density(m::MeshCartesianStatic, ACC0::Number)
    config = m.config
    accNewton = m.acc
    accNorm = sqrt.(accNewton.x .* accNewton.x + accNewton.y .* accNewton.y + accNewton.z .* accNewton.z)
    y = accNorm ./ ACC0
    RHS = similar(accNewton)
    RHS.x .= nu2.(y) .* accNewton.x
    RHS.y .= nu2.(y) .* accNewton.y
    RHS.z .= nu2.(y) .* accNewton.z

    #TODO MOND boundary conditions
    # Divergence of vector field
    ρx = diff_central_x(config.Δ[1], RHS.x)
    ρy = diff_central_y(config.Δ[2], RHS.y)
    ρz = diff_central_z(config.Δ[3], RHS.z)
    return ρPDM = ρx + ρy + ρz
end

function mesh_set_boundary_potential_QUMOND(index, m, A, b, TA, sq, ms)
    b[index] = sq * log10(norm(ustrip(m.pos[index]) - ms))
    A[index,:] .= zero(TA)
    A[index,index] = one(TA)
end

"""
    QUMOND_phi(m::MeshCartesianStatic, ACC0::Number, G::Number)

First compute ρ_PDM, then solve QUMOND (QUasi-linear MOdified Newtonian Dynamics) equation on the mesh.
Return modified potential
"""
function QUMOND_phi(m::MeshCartesianStatic, ACC0::Number, G::Number, Device::DeviceType, sparse::Bool)
    config = m.config

    ms = ustrip(mass_center(m))
    M = ustrip(total_mass(m))
    sq = ustrip(sqrt(G * M * ACC0))

    g = 4 * pi * G
    ρPDM = QUMOND_PDM_density(m, ACC0) .* g

    NX = config.N[1] + 2 * config.NG + 1
    NY = config.N[2] + 2 * config.NG + 1
    NZ = config.N[3] + 2 * config.NG + 1

    A = delta_mat3(NX, NY, NZ, ustrip.(getuLength(config.units), config.Δ)...;
        m.config.boundary, sparse,
    )
    b = ustrip.(Array(ρPDM .* (4 * pi * G))[:])

    # Manually set boundary potential and set diagonal element to one
    # Because the mesh can have no particle data, we compute potential by mesh.rho
    TA = eltype(A)
    Tb = eltype(b)
    Tphi = eltype(m.phi)

    # Omit edge points
    mesh_set_boundary_potential_QUMOND(1, m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX, m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX*(NY-1)+1, m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX*NY, m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(1 + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX*(NY-1)+1 + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
    mesh_set_boundary_potential_QUMOND(NX*NY*NZ, m, A, b, TA, sq, ms)

    for i in 2:NX-1
        mesh_set_boundary_potential_QUMOND(i, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*(NY-1)+i, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(i + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*(NY-1)+i + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
    end

    for j in 2:NY-1
        mesh_set_boundary_potential_QUMOND(NX*(j-1)+1, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*j, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*(j-1)+1 + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*j + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
    end

    for k in 2:NZ-1
        mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + 1, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + NX, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + NX*(NY-1)+1, m, A, b, TA, sq, ms)
        mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + NX*NY, m, A, b, TA, sq, ms)
    end

    # face
    for j in 2:NY-1
        for i in 2:NX-1
            mesh_set_boundary_potential_QUMOND(NY*(j-1)+i, m, A, b, TA, sq, ms)
            mesh_set_boundary_potential_QUMOND(NY*(j-1)+i + NX*NY*(NZ-1), m, A, b, TA, sq, ms)
        end
    end

    for k in 2:NZ-1
        for i in 2:NX-1
            mesh_set_boundary_potential_QUMOND(NX*NY*(k-1)+i, m, A, b, TA, sq, ms)
            mesh_set_boundary_potential_QUMOND(NX*NY*(k-1)+i + NX*(NY-1), m, A, b, TA, sq, ms)
        end
    end

    for k in 2:NZ-1
        for j in 2:NY-1
            mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + (j-1)*NX + 1, m, A, b, TA, sq, ms)
            mesh_set_boundary_potential_QUMOND(NX*NY*(k-1) + (j-1)*NX + NX, m, A, b, TA, sq, ms)
        end
    end

    return reshape(solve_matrix_equation(A, b, Device), NX, NY, NZ) .* unit(eltype(m.phi))
end

"""
    QUMOND_acc(m::MeshCartesianStatic, ACC0::Number, G::Number)

1. Compute ρ_PDM
2. Solve modified potential on the mesh
3. Compute acceleration by finite differencing the potential

Return acceleration
"""
function QUMOND_acc(m::MeshCartesianStatic, ACC0::Number, G::Number, Device::DeviceType, sparse::Bool)
    config = m.config
    phi_PDM = QUMOND_phi(m, ACC0, G, Device, sparse)
    acc = similar(m.acc)
    acc.x .= -diff_central_x(config.Δ[1], phi_PDM)
    acc.y .= -diff_central_y(config.Δ[2], phi_PDM)
    acc.z .= -diff_central_z(config.Δ[3], phi_PDM)
    return acc
end

"""
    QUMOND(m::MeshCartesianStatic, ACC0::Number, G::Number)

Apply QUMOND (QUasi-linear MOdified Newtonian Dynamics) formula to accelerations
"""
function QUMOND_acc!(m::MeshCartesianStatic, ACC0::Number, G::Number, Device::DeviceType, sparse::Bool)
    config = m.config
    phi_PDM = QUMOND_phi(m, ACC0, G, Device, sparse)
    m.acc.x .-= diff_central_x(config.Δ[1], phi_PDM)
    m.acc.y .-= diff_central_y(config.Δ[2], phi_PDM)
    m.acc.z .-= diff_central_z(config.Δ[3], phi_PDM)
    return nothing
end