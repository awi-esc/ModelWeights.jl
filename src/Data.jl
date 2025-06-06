module Data

export LEVEL, MODEL_MEMBER_DELIM
export DataMap, ClimateData

using CSV
using DataFrames
using Dates
using DimensionalData
using GLM
using NCDatasets
using JLD2
using Serialization
using Setfield
using Statistics
using YAML
using YAXArrays

# using GLM
# using Setfield

const MODEL_MEMBER_DELIM = "#"
const MODEL_NAME_FIXES = Dict("FGOALS_g2" => "FGOALS-g2", "ACCESS1.3" => "ACCESS1-3")


@enum LEVEL MODEL = 0 MEMBER = 1

KEYS_METADATA = [
    "mip_era",
    "variant_label",
    "grid_label",
    "realization",
    "physics_version",
    "initialization_method",
    "units",
]

const DataMap = Dict{String, YAXArray}

function Base.show(io::IO, x::Dict{String,YAXArray})
    println(io, "$(typeof(x))")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::Dict{String,YAXArray})
    println(io, "$(typeof(x))")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end


@kwdef mutable struct ClimateData
    info::Dict
    models::DataMap
    obs::DataMap
    weights::Union{YAXArray, Nothing}
end


function ClimateData(models::DataMap, observations::DataMap)
    return ClimateData(info=Dict(), models=models, obs=observations, weights=nothing)
end


function Base.show(io::IO, ::MIME"text/plain", x::ClimateData)
    println(io, "$(typeof(x))")
    models = collect(keys(x.models))
    obs = collect(keys(x.obs))
    shared = intersect(models, obs)    
    println("Model AND observational data: $(length(shared))")
    map(println, shared)
    println("----")
    println("Data only for observations: $(length(filter(x-> !(x in models), obs)))")
    println("Data only for models: $(length(filter(x-> !(x in obs), models)))")
    println("----")
    if !isnothing(x.weights)
        println("Weights: $(x.weights)")
    end
end

include("data-utils.jl")
include("data-functions.jl")
include("datamap.jl")
include("helper-functions.jl")
include("diagnostics.jl")

end