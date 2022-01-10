function apply_background_potential(f::Function, data::StructArray)
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds data.Potential[i] += f(data.Pos[i])
        end
    end
end