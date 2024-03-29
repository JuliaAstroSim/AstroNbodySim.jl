function cnn_poisson(m, u, dataML)
    unit_rho = m.rho / u.M * u.L^3
    unit_phi = dataML.dev_cpu(Lux.apply(
        dataML.tstate.model,
        dataML.dev_gpu(reshape(unit_rho, size(unit_rho)..., 1, 1)),
        dataML.tstate.parameters,
        dataML.tstate.states
    )[1][:,:,:,1,1])
    m.phi .= unit_phi * u.unitL^2 / unitT^2
end