struct GravToEvaluate{T<:AbstractPoint, S<:Number, L<:Number}
    Index::Int64
    Pos::T
    Collection::Collection
    OldAccTol::S
    SoftLen::L
end

struct GravResult{T<:AbstractPoint}
    Index::Int64
    Acc::T
    NInteractions::Int64
end

# tree potential uses GravToEvaluate, because OldAcc is vital for tree node opening
struct PotentialToEvaluate{T<:AbstractPoint, L<:Number}
    Index::Int64
    Pos::T
    Collection::Collection
    SoftLen::L
end

struct PotentialResult{PPM<:Number}
    Index::Int64
    Potential::PPM
end

function extract_grav(grav::GravToEvaluate)
    return grav.Pos, grav.Collection, grav.OldAccTol
end

function pack_grav_result(grav::GravToEvaluate, acc::AbstractPoint, ninteractions::Int64)
    return GravResult(grav.Index, acc, ninteractions)
end

function pack_pot_result(grav::GravToEvaluate, pot::Number)
    return PotentialResult(grav.Index, pot)
end

mutable struct Buffer{POS, ACC, Len, ACC0, PPM}
    #! For better performance, do not use buffer of abstract type
    sendbuffer::Dict{Int64, Any}
    recvbuffer::Dict{Int64, Any}

    gravToEval::Dict{Int64, Vector{GravToEvaluate{POS, ACC0, Len}}}
    gravToEvalRecv::Dict{Int64, Vector{GravToEvaluate{POS, ACC0, Len}}}
    gravResult::Dict{Int64, Vector{GravResult{ACC}}}
    gravResultRecv::Dict{Int64, Vector{GravResult{ACC}}}

    potToEval::Dict{Int64, Vector{PotentialToEvaluate{POS, Len}}}
    potToEvalRecv::Dict{Int64, Vector{PotentialToEvaluate{POS, Len}}}
    potResult::Dict{Int64, Vector{PotentialResult{PPM}}}
    potResultRecv::Dict{Int64, Vector{PotentialResult{PPM}}}
end

Buffer(PosType, AccType, LenType, OldAccType, PotType) = Buffer(
    Dict{Int64, Any}(),
    Dict{Int64, Any}(),
    Dict{Int64, Vector{GravToEvaluate{PosType, OldAccType, LenType}}}(),
    Dict{Int64, Vector{GravToEvaluate{PosType, OldAccType, LenType}}}(),
    Dict{Int64, Vector{GravResult{AccType}}}(),
    Dict{Int64, Vector{GravResult{AccType}}}(),

    Dict{Int64, Vector{PotentialToEvaluate{PosType, LenType}}}(),
    Dict{Int64, Vector{PotentialToEvaluate{PosType, LenType}}}(),
    Dict{Int64, Vector{PotentialResult{PotType}}}(),
    Dict{Int64, Vector{PotentialResult{PotType}}}(),
)

function Buffer(ZeroValues::ZeroValue)
    PosType = typeof(ZeroValues.pos)
    AccType = typeof(ZeroValues.acc)
    LenType = typeof(ZeroValues.pos.x)
    OldAccType = typeof(ZeroValues.acc.x)
    PotType = typeof(ZeroValues.potpermass)
    return Buffer(PosType, AccType, LenType, OldAccType, PotType)
end

Buffer(config::SimConfig) = Buffer(config.ZeroValues)