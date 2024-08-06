"""
$TYPEDSIGNATURES

Compute ρ_PDM on the RHS (right hand side) of QUMOND (QUasi-linear MOdified Newtonian Dynamics).
Return ρ_PDM
"""
function QUMOND_PDM_density(m::MeshCartesianStatic, ACC0::Number, G::Number)
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
    return ρPDM = (ρx + ρy + ρz) ./ (4π*G)
end

function mesh_set_boundary_potential_QUMOND(index, m, A, b, TA, sq, ms)
    b[index] = sq * log10(norm(ustrip(m.pos[index]) - ms))
    A[index,:] .= zero(TA)
    A[index,index] = one(TA)
end

"""
$TYPEDSIGNATURES

First compute ρ_PDM, then solve QUMOND (QUasi-linear MOdified Newtonian Dynamics) equation on the mesh.
Return modified potential

Warning: Currently the Poisson equation is solved with periodic boundary conditions
"""
function QUMOND_phi(m::MeshCartesianStatic, ACC0::Number, G::Number)
    rho_PDM = QUMOND_PDM_density(m, ACC0, G)
    rho = rho_PDM + m.rho
    fft_poisson(m.config.Δ, m.config.Len, ustrip.(rho*4π*G), Periodic(), m.config.device) .* unit(eltype(m.phi))
end

"""
$TYPEDSIGNATURES

1. Compute ρ_PDM
2. Solve modified potential on the mesh
3. Compute acceleration by finite differencing the potential

Return acceleration
"""
function QUMOND_acc(m::MeshCartesianStatic, ACC0::Number, G::Number)
    config = m.config
    phi_PDM = QUMOND_phi(m, ACC0, G)
    acc = similar(m.acc)
    acc.x .= -diff_central_x(config.Δ[1], phi_PDM)
    acc.y .= -diff_central_y(config.Δ[2], phi_PDM)
    acc.z .= -diff_central_z(config.Δ[3], phi_PDM)
    return acc
end

"""
$TYPEDSIGNATURES

Apply QUMOND (QUasi-linear MOdified Newtonian Dynamics) formula to accelerations
"""
function QUMOND_acc!(m::MeshCartesianStatic, ACC0::Number, G::Number)
    config = m.config
    phi_PDM = QUMOND_phi(m, ACC0, G)
    m.acc.x .-= diff_central_x(config.Δ[1], phi_PDM)
    m.acc.y .-= diff_central_y(config.Δ[2], phi_PDM)
    m.acc.z .-= diff_central_z(config.Δ[3], phi_PDM)
    return nothing
end