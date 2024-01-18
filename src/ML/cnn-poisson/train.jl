function train_cnn_poisson2d(;
    path = "output",
    epochs = 50, # For each dataset, train for `epoch` times, and step to next dataset
    train_nums = 10, # times that the train process will iterate
    lr = 0.001, # learning rate
    image_size = 16,
    data_nums = 100000, # training data
    batch_size = 50,
    read_data = false, # if false, generate new data for each run; if true, train with data in the last train
    read_model = false, # if true, continue train model based on the last train
    steps_between_plot = 50, # When plotting, the model will also be saved
    Realtime = false,
    resolution = (1600, 900)
)
    mkpathIfNotExist(path)
    mkpathIfNotExist(joinpath(path , "MSEloss"))
    mkpathIfNotExist(joinpath(path , "Maxloss"))
    mkpathIfNotExist(joinpath(path , "Result"))

    if isfile(joinpath(path, "stop"))
        rm(joinpath(path, "stop"), force=true)
    end

    if read_model==true
        @info("load model...")
        model = Knet.load(joinpath(path , "trained_model.jld2"), "model")
    else
        @info("build new model...")
        model = Chain((
            Conv2D(3,1,28),
            IdentitySkip2D(3,28,Conv2D),
            Conv2D(3,28,64),
            IdentitySkip2D(3,64,Conv2D),
            Conv2D(3,64,128),
            IdentitySkip2D(3,128,Conv2D),
            Conv2D(3,128,256),
            IdentitySkip2D(3,256,Conv2D),
            IdentitySkip2D(3,256,deConv2D),
            deConv2D(3,256,128),
            IdentitySkip2D(3,128,deConv2D),
            deConv2D(3,128,64),
            IdentitySkip2D(3,64,deConv2D),
            deConv2D(3,64,28),
            IdentitySkip2D(3,28,deConv2D),
            deConv2D(3,28,1),
        ))
    end

    if read_data==true && isfile(joinpath(path , "data_rho.jld2"))
        @info("Load data...")
        data_x = Knet.load(joinpath(path , "data_rho.jld2"), "data_x")
        data_y = Knet.load(joinpath(path , "data_phi.jld2"), "data_y")
    else
        @info("Training data inexist. Generate data...")
        data_x,data_y = Dataset(data_nums,image_size,batch_size)
        Knet.save(joinpath(path , "data_rho.jld2"), "data_x", data_x)
        Knet.save(joinpath(path , "data_phi.jld2"), "data_y", data_y)
    end
    
    MSElosses = Real[]
    Maxlosses = Real[]
    steps = Int[]
    
    if Realtime
        @info("Start visualization...")
        fig = Figure(resolution = (1600, 900))
        AxisMSEloss = GLMakie.Axis(fig[1,1], title = "MSE Loss", xlabel = "iter")
        AxisMaxloss = GLMakie.Axis(fig[1,2], title = "Max loss", xlabel = "iter")
        
        MSElossNode = Observable(MSElosses)
        MaxlossNode = Observable(Maxlosses)
        stepsNode = Observable(steps)
        
        lines!(AxisMSEloss, stepsNode, MSElossNode)
        lines!(AxisMaxloss, stepsNode, MaxlossNode)
        display(fig)
    end

    # log loss
    file_loss = open(joinpath(path, "log_loss.csv"), "w")
    write(file_loss, "iter,MSEloss,Maxloss\n")

    
    @info("Start Training...")
    len_data = length(data_x)
    iter = 0
    prog = Progress(train_nums*len_data*epochs; color = :light_blue)
    try
        for n = 1:train_nums
            for j = 1:len_data
                for i=1:epochs
                    # manual interuption
                    if isfile(joinpath(path, "stop"))
                        break
                    end

                    x = data_x[j]
                    y = data_y[j]
                    adam!(model,[(x,y)]; lr)
                    MSEloss = model(x,y)
                    Maxloss = model(x,y,abs)

                    
                    if ((j-1)*epochs+i)%steps_between_plot == 0 || iter == 1
                        # visualize
                        if Realtime
                            push!(steps, iter)
                            push!(MSElosses, MSEloss)
                            push!(Maxlosses, Maxloss)

                            stepsNode[] = steps
                            MSElossNode[] = MSElosses
                            MaxlossNode[] = Maxlosses

                            xlims!(AxisMSEloss, (0, iter))
                            xlims!(AxisMaxloss, (0, iter))

                            if iter == 1
                                ylims!(AxisMSEloss, (0.0, MSEloss))
                                ylims!(AxisMaxloss, (0.0, Maxloss))
                            end
                        end

                        # save model
                        _savemodel(model, joinpath(path , "model"))
                    end
                    #=
                    _showimage(model(reshape(x[:,:,:,1],(image_size,image_size,1,1))), x[:,:,:,1], y[:,:,:,1], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "1")
                    _showimage(model(reshape(x[:,:,:,2],(image_size,image_size,1,1))), x[:,:,:,2], y[:,:,:,2], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "2")
                    _showimage(model(reshape(x[:,:,:,3],(image_size,image_size,1,1))), x[:,:,:,3], y[:,:,:,3], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "3")
                    _showimage(model(reshape(x[:,:,:,4],(image_size,image_size,1,1))), x[:,:,:,4], y[:,:,:,4], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "4")
                    _showimage(model(reshape(x[:,:,:,5],(image_size,image_size,1,1))), x[:,:,:,5], y[:,:,:,5], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "5")
                    _showimage(model(reshape(x[:,:,:,6],(image_size,image_size,1,1))), x[:,:,:,6], y[:,:,:,6], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "6")
                    _showimage(model(reshape(x[:,:,:,7],(image_size,image_size,1,1))), x[:,:,:,7], y[:,:,:,7], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "7")
                    _showimage(model(reshape(x[:,:,:,8],(image_size,image_size,1,1))), x[:,:,:,8], y[:,:,:,8], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "8")
                    _showimage(model(reshape(x[:,:,:,9],(image_size,image_size,1,1))), x[:,:,:,9], y[:,:,:,9], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "9")
                    _showimage(model(reshape(x[:,:,:,10],(image_size,image_size,1,1))), x[:,:,:,10], y[:,:,:,10], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "10")
                    _showimage(model(reshape(x[:,:,:,11],(image_size,image_size,1,1))), x[:,:,:,11], y[:,:,:,11], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "11")
                    _showimage(model(reshape(x[:,:,:,12],(image_size,image_size,1,1))), x[:,:,:,12], y[:,:,:,12], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "12")
                    _showimage(model(reshape(x[:,:,:,13],(image_size,image_size,1,1))), x[:,:,:,13], y[:,:,:,13], ((j-1)*epochs+i)%steps_between_plot, image_size, date3, "13")
                    =#
                    iter += 1
                    ProgressMeter.update!(prog, iter; showvalues = [
                        ("iter", iter),
                        ("train_num", n),
                        ("MSEloss", MSEloss),
                        ("Maxloss", Maxloss),
                    ])
                    write(file_loss, @sprintf("%d,%.12f,%.12f\n", iter, MSEloss, Maxloss))
                end
            end
        end
    catch e
        throw(e)
        close(file_loss)
    end
    close(file_loss)

    if Realtime
        xlims!(AxisMSEloss, (0, iter))
        xlims!(AxisMaxloss, (0, iter))
        Makie.save(joinpath(path, "Train.png"), fig)
    else

    end
    @info("Training completed.")
end

function train_cnn_poisson3d(;
        
    )
    
end