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

    # Divergence of vector field
    # ρx = diff_central_x(config.Δ[1], RHS.x) 
    # ρy = diff_central_y(config.Δ[2], RHS.y) 
    # ρz = diff_central_z(config.Δ[3], RHS.z) 

    # Divergence of the ``nu2*acc`` vector field. ``PhysicalFDM.diff_central_*``
    # pads with zero at the boundary which, for a periodic mesh, is
    # inconsistent with the actual wrap-around and contaminates the
    # resulting ``rho_PDM`` whenever the RHS is not exactly zero at the
    # boundary. We instead compute the divergence with ``circshift``
    # which honors the periodic topology of ``m.config.boundary``.
    ρx = _periodic_central_diff(RHS.x, config.Δ[1], 1, config.boundary)
    ρy = _periodic_central_diff(RHS.y, config.Δ[2], 2, config.boundary)
    ρz = _periodic_central_diff(RHS.z, config.Δ[3], 3, config.boundary)
    return ρPDM = (ρx + ρy + ρz) ./ (4π*G)
end

"""
    _periodic_central_diff(u, Δ, axis, boundary)

Central finite difference of ``u`` along ``axis`` with spacing ``Δ``.
When ``boundary isa Periodic`` the wrap-around is honored via
``circshift`` (so a constant ``u`` has zero divergence everywhere);
otherwise fall back to ``diff_central_*`` which pads with zero.
"""
function _periodic_central_diff(u::AbstractArray, Δ, axis::Int, boundary::Periodic)
    shift_fwd = ntuple(i -> i == axis ? -1 : 0, ndims(u))
    shift_bwd = ntuple(i -> i == axis ? +1 : 0, ndims(u))
    return (circshift(u, shift_fwd) .- circshift(u, shift_bwd)) ./ (2 * Δ)
end

function _periodic_central_diff(u::AbstractArray, Δ, axis::Int, boundary)
    # Non-periodic: fall back to the standard ``diff_central_*`` from
    # PhysicalFDM which pads the boundary with zero.
    if axis == 1
        return diff_central_x(Δ, u)
    elseif axis == 2
        return diff_central_y(Δ, u)
    elseif axis == 3
        return diff_central_z(Δ, u)
    else
        error("Unsupported FDM dimension: $(axis)")
    end
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
    # ``m.rho`` is an ``ArrayScalarField``; convert to its underlying
    # data array so the broadcast ``rho_PDM + m.rho`` is plain
    # ``Array + Array``. This also avoids any ``RHOFIELD`` overload
    # ambiguity when ``m.rho`` shares its type parameter with
    # ``m.rho_e`` (PhysicalMeshes binds the two to one ``RHOFIELD``).
    rho = rho_PDM .+ m.rho.data
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
    acc.x .= .- _periodic_central_diff(phi_PDM, config.Δ[1], 1, config.boundary)
    acc.y .= .- _periodic_central_diff(phi_PDM, config.Δ[2], 2, config.boundary)
    acc.z .= .- _periodic_central_diff(phi_PDM, config.Δ[3], 3, config.boundary)
    return acc
end

"""
$TYPEDSIGNATURES

Apply QUMOND (QUasi-linear MOdified Newtonian Dynamics) formula to accelerations
"""
function QUMOND_acc!(m::MeshCartesianStatic, ACC0::Number, G::Number)
    config = m.config
    phi_PDM = QUMOND_phi(m, ACC0, G)
    m.acc.x .-= _periodic_central_diff(phi_PDM, config.Δ[1], 1, config.boundary)
    m.acc.y .-= _periodic_central_diff(phi_PDM, config.Δ[2], 2, config.boundary)
    m.acc.z .-= _periodic_central_diff(phi_PDM, config.Δ[3], 3, config.boundary)
    return nothing
end