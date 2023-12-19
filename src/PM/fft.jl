function fft_grid_kk(N, eps = 1e-6)
    hx = 2π / (N + 1)
    xx = Array{Float64}(undef, N + 1)
    
    # xx
    N_2 = floor(Int, 0.5 * N)
    for i in 1:N_2
        xx[i +  1] = hx * i
    end
    
    if isodd(N)
        N_2 = floor(Int, 0.5*(N + 1))
        for i in 1:N_2
            @inbounds xx[i + N_2] = hx * (i - N_2 - 1)
        end
    else
        N_2 = floor(Int, 0.5*(N + 1))
        for i in 1:N_2
            @inbounds xx[i + N_2 + 1] = hx * (i - N_2 - 1)
        end
    end
    
    xx[1] = eps # make sure that the denominator is not zero
    return xx
end

function fft_poisson(config::MeshConfig, rho::AbstractArray{T,1}, boundary::Periodic) where T
    rho_bar = fft(rho)    
    rho_bar[1] *= 0.0
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    xx = fft_grid_kk(config.Len[1])
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for i in 1:config.Len[1]+1
        u_bar[i] = rho_bar[i] / (delta2sum + delta2[1] * cos(xx[i]))
    end
    
    u = real(ifft(u_bar))
end

function fft_poisson(config::MeshConfig, rho::AbstractArray{T,2}, boundary::Periodic) where T
    rho_bar = fft(rho)    
    rho_bar[1] *= 0.0
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    xx = fft_grid_kk(config.Len[1])
    yy = fft_grid_kk(config.Len[2])
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for j in 1:config.Len[2]+1
        for i in 1:config.Len[1]+1
            u_bar[i,j] = rho_bar[i,j] / (delta2sum + delta2[1] * cos(xx[i]) + delta2[2] * cos(yy[j]))
        end
    end
    
    u = real(ifft(u_bar))
end

function fft_poisson(config::MeshConfig, rho::AbstractArray{T,3}, boundary::Periodic) where T
    rho_bar = fft(rho)    
    rho_bar[1] *= 0.0
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    xx = fft_grid_kk(config.Len[1])
    yy = fft_grid_kk(config.Len[2])
    zz = fft_grid_kk(config.Len[3])
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for k in 1:config.Len[3]+1
        for j in 1:config.Len[2]+1
            for i in 1:config.Len[1]+1
                u_bar[i,j,k] = rho_bar[i,j,k] / (delta2sum + delta2[1] * cos(xx[i]) + delta2[2] * cos(yy[j]) + delta2[3] * cos(zz[k]))
            end
        end
    end
    
    u = real(ifft(u_bar))
end

# Fast sine transform - homogeneous Dirichlet
function fft_poisson(config::MeshConfig, rho::AbstractArray{T,1}, boundary::Dirichlet) where T
    #rho_bar = fft(mesh.rho)
    #rho_bar[1] *= 0.0
    rho_bar = FFTW.r2r(complex(rho[2:end]), FFTW.RODFT00)
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    hx = pi / (config.Len[1] + 1)
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for i in 1:config.Len[1]
        u_bar[i] = rho_bar[i] / (delta2sum + delta2[1] * cos(hx * i))
    end

    u = real(FFTW.r2r(u_bar, FFTW.RODFT00)/((2*(config.Len[1] + 1))))
    #mesh.phi .= real(ifft(u_bar))
end

function fft_poisson(config::MeshConfig, rho::AbstractArray{T,2}, boundary::Dirichlet) where T
    #rho_bar = fft(mesh.rho)
    #rho_bar[1] *= 0.0
    rho_bar = FFTW.r2r(complex(rho[2:end, 2:end]), FFTW.RODFT00)
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    hx = pi / (config.Len[1] + 1)
    hy = pi / (config.Len[2] + 1)
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for j in 1:config.Len[2]
        for i in 1:config.Len[1]
            u_bar[i,j] = rho_bar[i,j] / (delta2sum + delta2[1] * cos(hx * i) + delta2[2] * cos(hy * j))
        end
    end

    u = real(FFTW.r2r(u_bar, FFTW.RODFT00)/((2*(config.Len[1] + 1)) * (2*(config.Len[2] + 1))))
    #mesh.phi .= real(ifft(u_bar))
end

function fft_poisson(config::MeshConfig, rho::AbstractArray{T,3}, boundary::Dirichlet) where T
    #rho_bar = fft(mesh.rho)
    #rho_bar[1] *= 0.0
    rho_bar = FFTW.r2r(complex(rho[2:end, 2:end, 2:end]), FFTW.RODFT00)
    
    delta2 = 2 ./ (ustrip.(config.Δ) .^ 2)
    delta2sum = - sum(delta2)
    
    hx = pi / (config.Len[1] + 1)
    hy = pi / (config.Len[2] + 1)
    hz = pi / (config.Len[3] + 1)
    
    # solve u_bar
    u_bar = similar(rho_bar);
    for k in 1:config.Len[3]
        for j in 1:config.Len[2]
            for i in 1:config.Len[1]
                u_bar[i,j,k] = rho_bar[i,j,k] / (delta2sum + delta2[1] * cos(hx * i) + delta2[2] * cos(hy * j) + delta2[3] * cos(hz * k))
            end
        end
    end

    u = real(FFTW.r2r(u_bar, FFTW.RODFT00)/((2*(config.Len[1] + 1)) * (2*(config.Len[2] + 1)) * (2*(config.Len[3] + 1))))
    #mesh.phi .= real(ifft(u_bar))
end

# Dirichlet BC returns a smaller array
function fft_poisson!(m::MeshCartesianStatic, rho::AbstractArray, boundary::Periodic)
    m.phi .= fft_poisson(m.config, rho, boundary) .* unit(eltype(m.phi))
end

function fft_poisson!(m::MeshCartesianStatic, rho::AbstractArray{T,1}, boundary::Dirichlet) where T
    m.phi[2:end] .= fft_poisson(m.config, rho, boundary) .* unit(eltype(m.phi))
end

function fft_poisson!(m::MeshCartesianStatic, rho::AbstractArray{T,2}, boundary::Dirichlet) where T
    m.phi[2:end,2:end] .= fft_poisson(m.config, rho, boundary) .* unit(eltype(m.phi))
end

function fft_poisson!(m::MeshCartesianStatic, rho::AbstractArray{T,3}, boundary::Dirichlet) where T
    m.phi[2:end,2:end,2:end] .= fft_poisson(m.config, rho, boundary) .* unit(eltype(m.phi))
end

function fft_poisson(m::AbstractMesh, G::Number)
    g = 4 * pi * G
    fft_poisson!(m, ustrip.(m.rho .* g), m.config.boundary)
end