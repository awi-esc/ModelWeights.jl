module Weights

import StatsBase: countmap

using DimensionalData
using Distributions
using Interpolations
using JLD2
#using LinearAlgebra
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


@kwdef struct ClimWIP
    performance_distances::DimArray
    independence_distances::DimArray
    Di::YAXArray # generalized distances each model wrt performance
    Sij::YAXArray # generalized distances between pairs of models
    w::YAXArray # contains all three types of normalized weights (wI, wP, combined)
    w_members::YAXArray # weights distributed evenly across resp. model members
    config::ConfigWeights # metadata
end


function Base.show(io::IO, x::ClimWIP)
    println(io, "::$(typeof(x)):")
    for m in dims(x.w, :model)
        println(io, "$m: $(round(x.w[model = At(m)].data[1], digits=3))")
    end
end

include("weights-utils.jl")

end