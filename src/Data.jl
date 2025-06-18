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


const MODEL_MEMBER_DELIM = "#"
const MODEL_NAME_FIXES = Dict("FGOALS_g2" => "FGOALS-g2", "ACCESS1.3" => "ACCESS1-3")


@enum LEVEL MODEL = 0 MEMBER = 1

const LEVEL_LOOKUP = Dict(
    "model" => MODEL,
    "member" => MEMBER
)

META_CMIP5 = ["physics_version", "realization", "initialization_method"]
META_CMIP6 = ["mip_era", "variant_label", "grid_label"]

mutable struct MetaData
    paths::Vector{String}
    variable::String
    experiment::String
    alias::String
    timerange::String
    subdir::String
    id::String
    is_model_data::Bool
    variable_long::Union{String, Nothing}
    statistic::Union{String, Nothing}
    esmvaltoolPreproc::Union{String, Nothing}

    function MetaData(
        paths::Vector{String},
        variable::String,
        experiment::String,
        alias::String;
        timerange::String = "",
        subdir::String = "",
        id::String = "",
        is_model_data::Bool = true,
        variable_long::Union{String, Nothing} = nothing,
        statistic::Union{String, Nothing} = nothing,
        esmvaltoolPreproc::Union{String, Nothing} = nothing
    )
        # @info "Constructor called."
        if any(map(isempty, [variable, experiment, alias]))
            throw(ArgumentError("alias (name of diagnostic in ESMValTool), variable and experiment must be provided!"))
        end
        absent(x) = isnothing(x) || isempty(x) 
        if absent(subdir)
            subdir = !absent(statistic) ? join([variable, statistic], "_") :
                throw(ArgumentError("Either subdir or statistic must be provided!"))
        end
        id_built = join([subdir, alias], "_")
        if !isempty(id)
            # if id and statistics are given, it must hold: id = variable_statistic
            if id_built != id
                throw(ArgumentError("ID must be identical to SUBDIR_ALIAS, if SUBDIR isn't provide it is assumed to equal VARIABLE_STATISTIC"))
            end
        end
        return new(
            paths,
            variable,
            experiment,
            alias,
            timerange,
            subdir,
            id_built,
            is_model_data,
            variable_long,
            statistic,
            esmvaltoolPreproc
        )
    end
end

# function Base.show(io::IO, x::MetaData)
#     println(io, "$(typeof(x)): $(x.id): ($(length(x.paths)) files)")
# end


const DataMap = Dict{String, YAXArray}
const MetaDataMap = Dict{String, MetaData}

function Base.show(io::IO, x::Dict{String, YAXArray})
    println(io, "$(typeof(x))")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::Dict{String, YAXArray})
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