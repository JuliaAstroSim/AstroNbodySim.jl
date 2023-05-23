# Generate differential matrices
#   Reference Credit: Qian, Long. 2021-09 (龙潜，github: longqian95)
#     - smooth_coef, diff_mat
#     - laplace_conv_op

"""
    smooth_coef(x_coord, fitting_order, diff_order)

Fit the data with `fitting_order` polynomial and calculate the `diff_order` differential at `x=0`. Input the x coordinates of data. Return the coefficients for calculating the result from data.

Suppose the data is `[u,v,w]`, we want to smooth `v` by fitting a line (1 order polynomial) and use it estimate a better `v`. Then set the x coordinates of `[u,v,w]` be `[-1,0,1]`, fit it and calculate the result (0 order differential) at x=0, the result will be `(u+v+w)/3`. So `smooth_coef([-1,0,1],1,0)` will return `[1/3,1/3/1/3]`.

Using `Rational` type or `Sym` type of `SymPy` can get the exact coefficients. For example, `smooth_coef(Sym[0,1,2],2,1)` get `Sym[-3/2, 2, -1/2]`, which means the first order differential of `[u,v,w]` at `u` is `-1.5u+2v-0.5w`. Using this way can generate all data in <https://en.wikipedia.org/wiki/Finite_difference_coefficient>

# Examples

- Linear extrapolation: `smooth_coef([1,2],1,0)` get `[2, -1]`, so the left linear extrapolation of `[u,v]` is `2u-v`.
- Quadratic extrapolation: `smooth_coef([-3,-2,-1],2,0)` get `[1, -3, 3]`, so the right Quadratic extrapolation of `[u,v,w]` is `u-3v+3w`.
- First order central differential: `smooth_coef(Rational[-1,0,1],2,1)` get `[-1//2, 0//1, 1//2]`, so the first order differential of `[u,v,w]` at v is `(w-u)/2`
- Second order central differential: `smooth_coef([-1,0,1],2,2)` get `[1, -2, 1]`, so the second order differential of `[u,v,w]` at v is `u+w-2v`
- Five points quadratic smoothing: `smooth_coef([-2,-1,0,1,2],2,0)` get `[-3/35, 12/35, 17/35, 12/35, -3/35]`, so `[-3/35 12/35 17/35 12/35 -3/35]*[a,b,c,d,e]` get the smoothed `c`.
- Four points quadratic interpolation: `smooth_coef([-3,-1,1,3],2,0)` get `[-1/16, 9/16, 9/16, -1/16]`, so `smooth_coef([-3,-1,1,3],2,0)'*[a,b,c,d]` get the estimated value at the middle of `b` and `c`.
"""
function smooth_coef(x_coord, fitting_order, diff_order)
    fitting_order >= length(x_coord) && throw(ArgumentError("cannot fit polynomial because length of $x_coord is not greater than the fitting order $fitting_order"))
    k=cat((x_coord.^i for i=0:fitting_order)...; dims=2)
    kt=transpose(k)
    p=inv(kt*k)*kt #Pseudo-inverse
    if diff_order<=fitting_order
        return factorial(diff_order)*p[diff_order+1,:]
    else
        return zeros(eltype(x_coord),length(x_coord))
    end
end

"""
    diff_mat(n, order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, lpoints=div(points,2), rpoints=points-lpoints-1, fitting_order=lpoints+rpoints, boundary=:Extrapolation, boundary_points=lpoints+rpoints+1, boundary_order=boundary_points-1, sparse=false)

Generate differential matrix.

# Arguments

- `n`: Matrix size. If `v` is length `n` vector, diff_mat(n)*v calculate the differential of `v`.
- `order`: Differential order. 
- `T`: Matrix element type. If set `T=Rational` or `using SymPy` and set `T=Sym`, `diff_mat` will return the exact value.
- `dt`: Numerical differential step size.
- `points`: Number of points for fitting polynomial to estimate differential. This argument is only for convenience. The real number of points is always `lpoints+rpoints+1`.
- `lpoints`: Number of points at left to the target point, which differential is calculated by the fitted polynomial.
- `rpoints`: Number of points at right to the target point. If `lpoints==rpoints`, then the differential is estimated as central finite difference. If `lpoints==0`, then it is normal forward finite difference. If `rpoints==0`, then it is backward finite difference.
- `fitting_order`: The order of the fitted polynomial for estimating differential.
- `boundary`: Boundary condition. Can be `Dirichlet()`(boundary value is zero), `Periodic`(assume data is periodic), `:Extrapolation`(boundary value is extrapolated according to `boundary_points` and `boundary_order`), `:None`(not deal with boundary, will return non-square matrix).
- `boundary_points`: Number of points for fitting polynomial to estimate differential at boundary. Normally it should not be much less than `points`, otherwise sometimes the current point may not be used to estimate the differential.
- `boundary_order`: The order of the fitted polynomial for points determined by `boundary_points`.
- `sparse`: If true, return sparse matrix instead of dense one.

# Examples

```
k=5; x=rand(k);
diff_mat(k,1;points=3)*x #do 3 points 1st order central differential ((x[n+1]-x[n-1])/2).
diff_mat(k,2;points=3)*x #do 3 points 2nd order central differential (x[n+1]+x[n-1]-2x[n]).
diff_mat(k,1;points=2,lpoints=0)*x #do the normal 1st order forward differential (x[n+1]-x[n]).
diff_mat(k,1;lpoints=1,rpoints=0)*x #do the 1st order backward differential (x[n-1]-x[n]).
```
"""
function diff_mat(n, order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, lpoints=div(points,2), rpoints=points-lpoints-1, fitting_order=lpoints+rpoints, boundary=Dirichlet(), boundary_points=lpoints+rpoints+1, boundary_order=boundary_points-1, sparse=false)
    n<lpoints+rpoints+1 && throw(ArgumentError("matrix size $n must be greater than or equal to lpoints+rpoints+1 ($lpoints+$rpoints+1)"))
    x=-lpoints:rpoints
    dt=dt^order
    v=smooth_coef(T.(x),fitting_order,order)/dt
    if sparse
        diagm1=spdiagm
        zeros1=spzeros
    else
        diagm1=diagm
        zeros1=zeros
    end
    if boundary isa Dirichlet || boundary isa Vacuum # in the vacuum case, we manually compute solution on the boundaries
        m=diagm1((x[i]=>v[i]*ones(T,n-abs(x[i])) for i=1:length(x))...)
    elseif boundary isa Periodic
        m=zeros1(T,n,n)
        for i=1:n, j=1:length(x)
            jj=mod1(x[j]+i,n)
            m[i,jj]=v[j]
        end
    #elseif boundary==:None #TODO
    #    nn=n-lpoints-rpoints
    #    m=diagm1(nn, n, (x[i]+lpoints=>v[i]*ones(T,nn) for i=1:length(x))...)
    #elseif boundary==:Extrapolation #TODO
    #    m=zeros1(T,n,n)
    #    for i=1:n
    #        if i<=lpoints
    #            b=smooth_coef(T.(1-i:boundary_points-i),boundary_order,order)/dt
    #            m[i,1:length(b)]=b
    #        elseif i>n-rpoints
    #            b=smooth_coef(T.(n-i-boundary_points+1:n-i),boundary_order,order)/dt
    #            m[i,n-length(b)+1:n]=b
    #        else
    #            m[i,x.+i]=v
    #        end
    #    end
    else
        throw(ArgumentError("unsupported boundary condition: $boundary"))
    end
    return m
end
diff_vec(order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, lpoints=div(points,2), rpoints=points-lpoints-1)=diff_mat(lpoints+rpoints+1,order;T=T,dt=dt,points=lpoints+rpoints+1,lpoints=lpoints,rpoints=rpoints,fitting_order=lpoints+rpoints,boundary=:None)[:]


#generate given order differential matrix for a vector which is expanded from row*col matrix
function diff_mat2_x(row,col,order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, boundary=Dirichlet(), sparse=false)
    t=diff_mat(col,order; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
    m=kron(t,I(row))
    return m
end
function diff_mat2_y(row,col,order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, boundary=Dirichlet(), sparse=false)
    t=diff_mat(row,order; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
    m=kron(I(col),t)
    return m
end

#2D Δ(Laplacian) operator
function delta_mat2(row,col; T=Float64, dt=one(T), points=3, boundary=Dirichlet(), sparse=false)
    return diff_mat2_x(row,col,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)+diff_mat2_y(row,col,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
end

function delta_mat2(row,col,Δx,Δy; T=Float64, dt=one(T), points=3, boundary=Dirichlet(), sparse=false)
    return diff_mat2_x(row,col,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)/Δx^2+diff_mat2_y(row,col,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)/Δy^2
end

#generate given order differential matrix for a vector which is expanded from row*col*page tensor
function diff_mat3_x(row,col,page,order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, boundary=Dirichlet(), sparse=false)
    t=diff_mat(col,order; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
    m=kron(I(page),kron(t,I(row)))
    return m
end
function diff_mat3_y(row,col,page,order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, boundary=Dirichlet(), sparse=false)
    t=diff_mat(row,order; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
    m=kron(I(col*page),t) #or: m=kron(I(page),kron(I(col),t))
    return m
end
function diff_mat3_z(row,col,page,order=1; T=Float64, dt=one(T), points=2*div(order+1,2)+1, boundary=Dirichlet(), sparse=false)
    t=diff_mat(page,order; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
    m=kron(t,I(row*col)) #or: m=kron(kron(t,I(row)),I(col))
    return m
end


#3D Δ(Laplacian) operator
function delta_mat3(row,col,page; T=Float64, dt=one(T), points=3, boundary=Dirichlet(), sparse=false)
    return diff_mat3_x(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)+diff_mat3_y(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)+diff_mat3_z(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)
end

function delta_mat3(row,col,page,Δx,Δy,Δz; T=Float64, dt=one(T), points=3, boundary=Dirichlet(), sparse=false)
    return diff_mat3_x(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)/Δx^2+diff_mat3_y(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)/Δy^2+diff_mat3_z(row,col,page,2; T=T,dt=dt,points=points,boundary=boundary,sparse=sparse)/Δz^2
end

function conv(kernel::AbstractArray{T,1}, d::AbstractArray{T,1}, boundary=Dirichlet(); fill = zero(T)) where T
    h = div(length(kernel), 2)
    d1 = PaddedView(fill, d, (1-h:length(d)+h,))
    @tullio out[x] := d1[x-i] * kernel[i]
    return out
end

function conv(kernel::AbstractArray{T,2}, d::AbstractArray{T,2}, boundary=Dirichlet(); fill = zero(T)) where T
    h=div.(size(kernel),2)
    d1=PaddedView(fill,d,(1-h[1]:size(d,1)+h[1],1-h[2]:size(d,2)+h[2]))
    @tullio out[x,y]:=d1[x-i,y-j]*kernel[i,j]
    return parent(out)
end

function conv(kernel::AbstractArray{T,3}, d::AbstractArray{T,3}, boundary=Dirichlet(); fill = zero(T)) where T
    h=div.(size(kernel),2)
    d1=PaddedView(fill,d,size(d).+h.+h,h.+1)
    @tullio out[x,y,z]:=d1[x-i,y-j,z-k]*kernel[i,j,k]
    return parent(out)
end

# ∇(partial) difference operator
function diff_oneside_op(Δx)
    SVector{2}([-1,1]) / Δx
end

function diff_central_op(Δx)
    SVector{3}([-1,0,1]) / Δx / 2.0
end

function diff2_central_op(Δx)
    SVector{3}([1,-2,1]) / Δx / Δx
end

function diff_central_x(Δx, u::AbstractArray{T,1}, pad = zero(T)) where T
    LenX = length(u)
    kernel = diff_central_op(Δx)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1-h:LenX+h,))
    @tullio out[x]:=d1[x+i] * kernel[i]
    return parent(out)
end

function diff_central_x(Δx, u::AbstractArray{T,2}, pad = zero(T)) where T
    LenX, LenY = size(u)
    kernel = diff_central_op(Δx)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1-h:LenX+h,     1:LenY))
    @tullio out[x,y]:=d1[x+i,y] * kernel[i]
    return parent(out)
end

function diff_central_y(Δy, u::AbstractArray{T,2}, pad = zero(T)) where T
    LenX, LenY, LenZ = size(u)
    kernel = diff_central_op(Δy)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1:LenX,     1-h:LenY+h))
    @tullio out[x,y]:=d1[x,y+i] * kernel[i]
    return parent(out)
end

function diff_central_x(Δx, u::AbstractArray{T,3}, pad = zero(T)) where T
    LenX, LenY, LenZ = size(u)
    kernel = diff_central_op(Δx)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1-h:LenX+h,     1:LenY, 1:LenZ))
    @tullio out[x,y,z]:=d1[x+i,y,z] * kernel[i]
    return parent(out)
end

function diff_central_y(Δy, u::AbstractArray{T,3}, pad = zero(T)) where T
    LenX, LenY, LenZ = size(u)
    kernel = diff_central_op(Δy)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1:LenX,     1-h:LenY+h, 1:LenZ))
    @tullio out[x,y,z]:=d1[x,y+i,z] * kernel[i]
    return parent(out)
end

function diff_central_z(Δz, u::AbstractArray{T,3}, pad = zero(T)) where T
    LenX, LenY, LenZ = size(u)
    kernel = diff_central_op(Δz)
    h = div(length(kernel), 2)
    d1 = PaddedView(pad, u, (1:LenX,     1:LenY, 1-h:LenZ+h))
    @tullio out[x,y,z]:=d1[x,y,z+i] * kernel[i]
    return parent(out)
end

# Δ(Laplacian) operator
laplace_conv_op(Δx) = diff2_central_op(Δx)

function laplace_conv_op(Δx, Δy)
    SMatrix{3,3}([0       1/Δy/Δy          0;
                  1/Δx/Δx -2/Δx/Δx-2/Δy/Δy 1/Δx/Δx;
                  0       1/Δy/Δy          0])
end

function laplace_conv_op(Δx, Δy, Δz)
    SArray{Tuple{3,3,3}}(cat(
    [0    0    0;
     0 1/Δz/Δz 0;
     0    0    0],
    [0       1/Δy/Δy                  0;
     1/Δx/Δx -2/Δx/Δx-2/Δy/Δy-2/Δz/Δz 1/Δx/Δx;
     0       1/Δy/Δy                  0],
    [0    0    0;
     0 1/Δz/Δz 0;
     0    0    0], dims = 3))
end

function laplace_conv(a, Δ...)
    kernel = laplace_conv_op(Δ...)
    return imfilter(a, kernel, Fill(0))
end

function solve_matrix_equation(A::Matrix, b, ::CPU)
    # TODO: use left devision
    # TODO: support units
    return pinv(A) * b
end

function solve_matrix_equation(A::SparseMatrixCSC, b, ::CPU)
    # TODO: support units
    return A \ b
end

function solve_matrix_equation(A, b, ::GPU)
    return cu(A) \ cu(b)
end

function fdm_poisson(m::AbstractMesh, dim::Val{1}, G::Number, ::VertexMode, Device::DeviceType, boundary::BoundaryCondition, sparse::Bool) # boundary: Periodic, Dirichlet
    config = m.config

    A = diff_mat(config.N[1] + 1 + 2 * config.NG, 2;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (config.Δ[1] * config.Δ[1] * 4 * pi * G)))
    
    phi = solve_matrix_equation(A, b, Device)
    m.phi .= phi * unit(eltype(m.phi))
end

function mesh_potential_1d(x, dx, rho, pos, G, Tphi)
    phi = zero(Tphi)
    for i in eachindex(pos)
        r = abs(pos[i] - x)
        if !iszero(r)
            phi -= dx * rho[i] * G / r
        end
    end
    return phi
end

function fdm_poisson(m::AbstractMesh, dim::Val{1}, G::Number, ::VertexMode, Device::DeviceType, boundary::Vacuum, sparse::Bool)
    config = m.config

    A = diff_mat(config.N[1] + 1 + 2 * config.NG, 2;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (config.Δ[1] * config.Δ[1] * 4 * pi * G)))
    
    # Manually set boundary potential (one cell outside the mesh)
    # Because the mesh can have no particle data, we compute potential by mesh.rho
    dx = config.Δ[1]
    b[1] -= mesh_potential_1d(m.pos[1] - dx, dx, m.rho, m.pos, G, eltype(m.phi))
    b[end] -= mesh_potential_1d(m.pos[end] + dx, dx, m.rho, m.pos, G, eltype(m.phi))
    
    phi = solve_matrix_equation(A, b, Device)
    m.phi .= phi * unit(eltype(m.phi))
end

function fdm_poisson(m::AbstractMesh, dim::Val{2}, G::Number, ::VertexMode, Device::DeviceType, boundary::BoundaryCondition, sparse::Bool) # boundary: Periodic, Dirichlet
    config = m.config
    A = delta_mat2((config.N .+ (1 + 2 * config.NG))..., ustrip.(getuLength(config.units), config.Δ)...;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (4 * pi * G)))[:]

    phi = solve_matrix_equation(A, b, Device)
    m.phi[:] .= phi * unit(eltype(m.phi))
end

function mesh_potential(p, Δ, rho, pos, G, Tphi)
    phi = zero(Tphi)
    #TODO multi-threading
    #TODO potential for different dimensions?
    for i in eachindex(pos)
        r = norm(pos[i] - p)
        if !iszero(r)
            @inbounds phi -= prod(Δ) * rho[i] * G / r
        end
    end
    return phi
end

function mesh_omit_corner(index, A, b, TA, Tb)
    A[index,:] .= zero(TA)
    A[index,index] = one(TA)
    b[index] = zero(Tb)
end

function mesh_set_boundary_potential(index, m, A, b, TA, Tphi, G)
    b[index] = ustrip(mesh_potential(m.pos[index], m.config.Δ, m.rho, m.pos, G, Tphi))
    A[index,:] .= zero(TA)
    A[index,index] = oneunit(TA)
end

function fdm_poisson(m::AbstractMesh, dim::Val{2}, G::Number, ::VertexMode, Device::DeviceType, boundary::Vacuum, sparse::Bool)
    config = m.config

    NX = config.N[1] + 2 * config.NG + 1
    NY = config.N[2] + 2 * config.NG + 1

    A = delta_mat2(NX, NY, ustrip.(getuLength(config.units), config.Δ)...;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (4 * pi * G)))[:]

    # Manually set boundary potential and set diagonal element to one
    # Because the mesh can have no particle data, we compute potential by mesh.rho
    TA = eltype(A)
    Tb = eltype(b)
    Tphi = eltype(m.phi)
    
    # cornor points
    mesh_set_boundary_potential(1, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*(NY-1)+1, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*NY, m, A, b, TA, Tphi, G)

    # edge
    for i in 2:NX-1
        mesh_set_boundary_potential(i, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*(NY-1)+i, m, A, b, TA, Tphi, G)
    end

    for j in 2:NY-1
        mesh_set_boundary_potential(NX*(j-1)+1, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*j, m, A, b, TA, Tphi, G)
    end

    phi = solve_matrix_equation(A, b, Device)
    m.phi[:] .= phi * unit(eltype(m.phi))
end

function fdm_poisson(m::AbstractMesh, dim::Val{3}, G::Number, ::VertexMode, Device::DeviceType, boundary::BoundaryCondition, sparse::Bool) # boundary: Periodic, Dirichlet
    config = m.config

    A = delta_mat3((config.N .+ (1 + 2 * config.NG))..., ustrip.(getuLength(config.units), config.Δ)...;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (4 * pi * G)))[:]

    phi = solve_matrix_equation(A, b, Device)
    m.phi[:] .= phi * unit(eltype(m.phi))
end

function fdm_poisson(m::AbstractMesh, dim::Val{3}, G::Number, ::VertexMode, Device::DeviceType, boundary::Vacuum, sparse::Bool)
    config = m.config

    NX = config.N[1] + 2 * config.NG + 1
    NY = config.N[2] + 2 * config.NG + 1
    NZ = config.N[3] + 2 * config.NG + 1

    A = delta_mat3(NX, NY, NZ, ustrip.(getuLength(config.units), config.Δ)...;
        boundary, sparse,
    )
    b = ustrip.(Array(m.rho .* (4 * pi * G)))[:]
    
    # Manually set boundary potential and set diagonal element to one
    # Because the mesh can have no particle data, we compute potential by mesh.rho
    TA = eltype(A)
    Tb = eltype(b)
    Tphi = eltype(m.phi)

    # Omit edge points
    mesh_set_boundary_potential(1, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*(NY-1)+1, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*NY, m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(1 + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*(NY-1)+1 + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
    mesh_set_boundary_potential(NX*NY*NZ, m, A, b, TA, Tphi, G)

    for i in 2:NX-1
        mesh_set_boundary_potential(i, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*(NY-1)+i, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(i + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*(NY-1)+i + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
    end

    for j in 2:NY-1
        mesh_set_boundary_potential(NX*(j-1)+1, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*j, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*(j-1)+1 + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*j + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
    end

    for k in 2:NZ-1
        mesh_set_boundary_potential(NX*NY*(k-1) + 1, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*NY*(k-1) + NX, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*NY*(k-1) + NX*(NY-1)+1, m, A, b, TA, Tphi, G)
        mesh_set_boundary_potential(NX*NY*(k-1) + NX*NY, m, A, b, TA, Tphi, G)
    end

    # face
    for j in 2:NY-1
        for i in 2:NX-1
            mesh_set_boundary_potential(NY*(j-1)+i, m, A, b, TA, Tphi, G)
            mesh_set_boundary_potential(NY*(j-1)+i + NX*NY*(NZ-1), m, A, b, TA, Tphi, G)
        end
    end

    for k in 2:NZ-1
        for i in 2:NX-1
            mesh_set_boundary_potential(NX*NY*(k-1)+i, m, A, b, TA, Tphi, G)
            mesh_set_boundary_potential(NX*NY*(k-1)+i + NX*(NY-1), m, A, b, TA, Tphi, G)
        end
    end

    for k in 2:NZ-1
        for j in 2:NY-1
            mesh_set_boundary_potential(NX*NY*(k-1) + (j-1)*NX + 1, m, A, b, TA, Tphi, G)
            mesh_set_boundary_potential(NX*NY*(k-1) + (j-1)*NX + NX, m, A, b, TA, Tphi, G)
        end
    end

    phi = solve_matrix_equation(A, b, Device)
    m.phi[:] .= phi * unit(eltype(m.phi))
end
