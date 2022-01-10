function substract_by_id(data::Array, symbol::Symbol; info::Symbol = Symbol())
    id = [p.ID for p in data]
    a = [getfield(p, symbol) for p in data]

    # substract additional info
    if sizeof(info) > 0
        ai = [getfield(p, info) for p in data]
        df = DataFrame(ID = id, value = a, ai = ai)
        rename!(df, [:ID, symbol, info])
    else
        df = DataFrame(ID = id, value = a)
        rename!(df, [:ID, symbol])
    end

    sort!(df)
    return df
end

function substract_by_id(data::StructArray, symbol::Symbol; info::Symbol = Symbol())
    # substract additional info
    if sizeof(info) > 0
        df = DataFrame(ID = data.ID, value = getproperty(data, symbol), ai = getproperty(data, info))
        rename!(df, [:ID, symbol, info])
    else
        df = DataFrame(ID = data.ID, value = getproperty(data, symbol))
        rename!(df, [:ID, symbol])
    end

    sort!(df)
    return df
end

function diff_by_id(df1::DataFrame, df2::DataFrame, symbol::Symbol; info::Symbol = Symbol())
    df = DataFrame(
        ID = empty(df1.ID),
        s1 = empty(getproperty(df1, symbol)),
        s2 = empty(getproperty(df1, symbol)),
        diff = empty(getproperty(df1, symbol)),
    )

    if sizeof(info) > 0
        df[!, info] = empty(getproperty(df1, info))
    end

    NumTotal = length(df1.ID)
    N = length(df2.ID)
    
    NumMatched = 0
    p = Progress(NumTotal)
    for i in 1:NumTotal
        id = df1[i, :ID]
        for k in 1:N
            if id == df2[k, :ID]
                s1 = df1[i, symbol]
                s2 = df2[k, symbol]
                if sizeof(info) > 0
                    push!(df, [id, s1, s2, s2 - s1, df1[i, info]])
                else
                    push!(df, [id, s1, s2, s2 - s1])
                end
                NumMatched += 1
                break
            end
        end
        next!(p)
    end

    if sizeof(info) > 0
        rename!(df, [:ID, Symbol(symbol, "_1"), Symbol(symbol, "_2"), :diff, info])
    else
        rename!(df, [:ID, Symbol(symbol, "_1"), Symbol(symbol, "_2"), :diff])
    end
    println("$NumTotal in total, $NumMatched ID matched")

    return df
end

function diff_by_id(data1, data2, symbol::Symbol; info::Symbol = Symbol())
    df1 = substract_by_id(data1, symbol; info)
    df2 = substract_by_id(data2, symbol; info)
    return diff_by_id(df1, df2, symbol; info)
end