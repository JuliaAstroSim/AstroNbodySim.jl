using Knet
using CUDA
using PhysicalFDM
using Printf
using GLMakie
using ProgressMeter
using Images


struct Chain
    layers
end
(c::Chain)(x) = (for l in c.layers;x = l(x); end; x)
(c::Chain)(x,y) = mean(abs2.(c(x).-y))
(c::Chain)(x,y,z) = maximum(z.(c(x).-y))



struct Conv2D
    w
    b
    a
end
(c::Conv2D)(x) = c.a .* relu.(conv4(c.w, x,padding=1)) .+ c.b
Conv2D(kernel_size::Int,in_channel::Int,out_channel::Int) = Conv2D(param(kernel_size,kernel_size,in_channel,out_channel;atype=Atype), param0(1,1,out_channel,1;atype=Atype), param(1,1,out_channel,1;atype=Atype))

struct deConv2D
    w
    b
    a
end
(d::deConv2D)(x) = d.a .* relu.(deconv4(d.w, x,padding=1)) .+ d.b
deConv2D(kernel_size::Int,in_channel::Int,out_channel::Int) = deConv2D(param(kernel_size,kernel_size,out_channel,in_channel;atype=Atype), param0(1,1,out_channel,1;atype=Atype), param(1,1,out_channel,1;atype=Atype))

struct pad0Conv2D
    w
    b
end
(c::pad0Conv2D)(x) = relu.(conv4(c.w, x,padding=0) .+ c.b)
pad0Conv2D(kernel_size,in_channel,out_channel) = pad0Conv2D(param(kernel_size,kernel_size,in_channel,out_channel;atype=Atype), param(1,1,out_channel,1;atype=Atype))

struct IdentitySkip2D
    inner
end
(m::IdentitySkip2D)(x) = m.inner(x) .+ x
IdentitySkip2D(kernel_size,keep_channel,f) = IdentitySkip2D(f(kernel_size,keep_channel,keep_channel))




struct Conv3D
    w
    b
    a
end
(c::Conv3D)(x) = c.a .* relu.(conv4(c.w, x,padding=1)) .+ c.b
Conv3D(kernel_size::Int,in_channel::Int,out_channel::Int) = Conv3D(param(kernel_size,kernel_size,kernel_size,in_channel,out_channel;atype=Atype), param0(1,1,1,out_channel,1;atype=Atype), param(1,1,1,out_channel,1;atype=Atype))

struct deConv3D
    w
    b
    a
end
(d::deConv3D)(x) = d.a .* relu.(deconv4(d.w, x,padding=1)) .+ d.b
deConv3D(kernel_size::Int,in_channel::Int,out_channel::Int) = deConv3D(param(kernel_size,kernel_size,kernel_size,out_channel,in_channel;atype=Atype), param0(1,1,1,out_channel,1;atype=Atype), param(1,1,1,out_channel,1;atype=Atype))

struct pad0Conv3D
    w
    b
end
(c::pad0Conv3D)(x) = relu.(conv4(c.w, x,padding=0) .+ c.b)
#TODO: correct?
pad0Conv3D(kernel_size,in_channel,out_channel) = pad0Conv3D(param(kernel_size,kernel_size,in_channel,out_channel;atype=Atype), param(1,1,out_channel,1;atype=Atype))

struct IdentitySkip3D
    inner
end
(m::IdentitySkip3D)(x) = m.inner(x) .+ x
IdentitySkip3D(kernel_size,keep_channel,f) = IdentitySkip3D(f(kernel_size,keep_channel,keep_channel))




function _savemodel(model, path)
    filename = joinpath(path, "trained_model.jld2")
    #@info("saving model to " * filename)
    Knet.save(filename, "model", model)
end
