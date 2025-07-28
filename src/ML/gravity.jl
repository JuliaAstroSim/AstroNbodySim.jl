function cnn_poisson(m, u, dataML)
    unit_rho = m.rho / u.M * u.L^3
    result = Lux.apply(
        dataML.tstate.model,
        dataML.dev_gpu(reshape(unit_rho, size(unit_rho)..., 1, 1)),
        dataML.tstate.parameters,
        dataML.tstate.states,
    )[1]
    unit_phi = reshape(unit_rho, size(unit_rho)...)
    m.phi .= unit_phi * u.L^2 / u.T^2
end