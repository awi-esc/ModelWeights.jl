module Weights

import StatsBase.ecdf

using Dates
using DimensionalData
using Distributions
using JLD2
using NCDatasets
using Serialization
using Setfield
using Statistics
using YAXArrays


using ..Data


@kwdef struct ConfigWeights
    performance::Dict{String, Number} = Dict()
    independence::Dict{String, Number} = Dict()
    sigma_performance::Number = 0.5
    sigma_independence::Number = 0.5
end

@kwdef struct MMEWeights
    w::YAXArray
    name::String

    function MMEWeights(w::YAXArray, name::String)
        Data.throwErrorIfDimMissing(w, :model)
        # the weight array must only have dimension :model
        if length(dims(w)) != 1
            throw(ArgumentError("MMEWEights array must only have dimension :model!"))
        end
        if isempty(name)
            @warn "MMEWeights wasn't assigned a name!"
        end
        return new(w, name)
    end
end

@kwdef struct ClimWIP
    performance_distances::Dict # from diagnostics to performance distances for all members
    independence_distances::Dict # from diagnostics to independence distances for all members
    Di::YAXArray # generalized distances each model wrt performance
    Sij::YAXArray # generalized distances between pairs of models
    w::YAXArray # contains all three types of normalized weights (wI, wP, combined)
    config::ConfigWeights # metadata
end


function Base.show(io::IO, x::ClimWIP)
    println(io, "::$(typeof(x)):")
    println("weights: $(Array(x.w.weight))")
    println("#models: $(length(x.w.model))")
end

include("weights-utils.jl")

end