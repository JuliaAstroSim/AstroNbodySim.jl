const Atype = KnetArray{Float32}
#const Atype = CuArray{Float32}

function make_data(data_nums,data_size,sorts)

    if (sorts==0)||(sorts==4)||(sorts==6)||(sorts==10)
        data = rand(Float64, data_size, data_size, 1, data_nums)
    end

    if sorts == 1
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[5:end-4,5:end-4,1,:] .= 0.0
    end

    if sorts == 3
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[1:end-10,11:end,1,:] .= 0.0
        data[11:end,1:end-10,1,:] .= 0.0
        data[1:end-10,1:end-10,1,:] .= 0.0
        data[11:end,11:end,1,:] .= 0.0
    end

    if sorts == 5
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[5:end-4,5:end-4,1,:] .= 0.1
    end

    if sorts == 7
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[1:end-10,11:end,1,:] .= 0.1
        data[11:end,1:end-10,1,:] .= 0.1
        data[1:end-10,1:end-10,1,:] .= 0.1
        data[11:end,11:end,1,:] .= 0.1
    end

    if sorts == 9
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[5:end-4,5:end-4,1,:] .= 0.99
    end

    if sorts == 11
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[1:end-10,11:end,1,:] .= 0.98
        data[11:end,1:end-10,1,:] .= 0.98
        data[1:end-10,1:end-10,1,:] .= 0.98
        data[11:end,11:end,1,:] .= 0.98
    end

    if (sorts==12)||(sorts==2)||(sorts==8)
        data = rand(Float64, data_size, data_size, 1, data_nums)
        data[4:end-3,4:end-3] = imfilter(data[4:end-3,4:end-3],Kernel.gaussian(2))
    end

    return data
end

function laplacian_DataGenerate(data_nums,data_size;atype=Atype)
    x = zeros(Float64, data_size, data_size, 1, data_nums)
    center = zeros(Float64,data_size-4,data_size-4,1,1)
    for i = 1:data_nums
        center = make_data(1,data_size-4,i%13)
        x[3:end-2,3:end-2,1,i] = center
    end
    y = zeros(Float64, data_size, data_size, 1, data_nums)
    for i = 1:data_nums
        grad2 = laplace_conv(x[:,:,1,i], 1.0, 1.0) #TODO add ??x, ??y
        y[:,:,1,i] = grad2
    end
    return convert.(atype,(y,x))
end


function make_minibatch(X, Y, idxs)
    X_batch = Atype(undef,size(X[:,:,:,1])..., length(idxs))
    for i in 1:length(idxs)
        X_batch[:,:,:,i] = X[:,:,:,idxs[i]]
    end
    Y_batch = Atype(undef,size(Y[:,:,:,1])..., length(idxs))
    for i in 1:length(idxs)
        Y_batch[:,:,:,i] = Y[:,:,:,idxs[i]]
    end
    return(X_batch, Y_batch)
end

function Dataset(data_nums,data_size,batch_size)
    X,Y = laplacian_DataGenerate(data_nums,data_size)
    mb_idxs = collect(partition(1:data_nums, batch_size))
    x,y=[],[]
    for i in mb_idxs
        append!(x, [make_minibatch(X, Y, i)[1]])
        append!(y, [make_minibatch(X, Y, i)[2]])
    end
    return (x,y)
end