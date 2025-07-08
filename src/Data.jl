module Data

export Level, MODEL_MEMBER_DELIM
export DataMap, ClimateData, ESMEnsemble

export subsetModelData, sharedModels, filterPathsSharedModels!
export summarizeEnsembleMembersVector, summarizeEnsembleMembersVector!
export alignPhysics, addMasks!
export defineDataMap, loadPreprocData


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
const MODEL_NAME_FIXES = Dict(
    "FGOALS_g2" => "FGOALS-g2", 
    "ACCESS1.3" => "ACCESS1-3",
    "fio-esm" => "FIO-ESM"
)

@enum Level MODEL_LEVEL = 0 MEMBER_LEVEL = 1
const LEVEL_LOOKUP = Dict(
    "model" => MODEL_LEVEL,
    "member" => MEMBER_LEVEL,
    :model => MODEL_LEVEL,
    :member => MEMBER_LEVEL
)
function getLevel(level::Union{String, Symbol})
    msg = "Model level must be one of $(keys(LEVEL_LOOKUP)), but found: $(level)."
    get(() -> throw(ArgumentError(msg)), LEVEL_LOOKUP, level)
    return level in [:model, "model"] ? MODEL_LEVEL : MEMBER_LEVEL
end


mutable struct MetaData
    variable::String
    experiment::String
    alias::String
    timerange::String
    subdir::String
    id::String
    variable_long::Union{String, Nothing}
    statistic::Union{String, Nothing}
    esmvaltoolPreproc::Union{String, Nothing}

    function MetaData(
        variable::String,
        experiment::String,
        alias::String,
        timerange::String,
        subdir::String,
        id::String,
        variable_long::Union{String, Nothing},
        statistic::Union{String, Nothing},
        esmvaltoolPreproc::Union{String, Nothing}
    )
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
            variable,
            experiment,
            alias,
            timerange,
            subdir,
            id_built,
            variable_long,
            statistic,
            esmvaltoolPreproc
        )
    end
end

function MetaData(
    variable::String,
    experiment::String,
    alias::String;
    timerange::String = "",
    subdir::String = "",
    id::String = "",
    variable_long::Union{String, Nothing} = nothing,
    statistic::Union{String, Nothing} = nothing,
    esmvaltoolPreproc::Union{String, Nothing} = nothing
)
    return MetaData(variable, experiment, alias, timerange, subdir, id, variable_long, statistic, esmvaltoolPreproc)
end


const DataMap = Dict{String, YAXArray}

function defineDataMap(data::Vector{<:YAXArray}, ids::Vector{String})
    if length(data) != length(ids)
        throw(ArgumentError("data and ids must have same size to build up DataMap"))
    end
    datamap = DataMap()
    for (d, id) in zip(data, ids)
        datamap[id] = d
    end
    return datamap
end


struct ClimateData
    models::DataMap
    obs::DataMap
end


const FILENAME_FORMATS = Dict(
    :esmvaltool => ("mip", "model", "table_id", "exp", "variant", "variable", "grid"),
    :cmip => ("variable", "table_id", "model",  "exp", "variant", "grid", "timerange")
)

@kwdef struct FilenameMeta
    variable::String
    table_id::String
    model::String
    experiment::String
    variant::String
    fn::String
    grid::Union{String, Nothing} = nothing # only necessary for CMIP6, not CMIP5
    timerange::Union{String, Nothing} = nothing # only given in filenames of CMIP standard
    mip::Union{String, Nothing} = nothing # in CMIP standard only given for CMIP5, not CMIP6
end


@kwdef mutable struct ESMEnsemble
    weights::Union{YAXArray, Nothing} 
    variables::Dict{Union{Symbol, String}, YAXArray}
    s::DataFrame   # styles for plotting
    p::DataFrame   # parameters
end


function Base.show(io::IO, x::Dict{String, YAXArray})
    println(io, "$(typeof(x))")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end

function Base.show(io::IO, x::ClimateData)
    println(io, "$(typeof(x))")
    println(io, "models: $(x.models)observations: $(x.obs)")
end

function Base.show(io::IO, x::Union{MetaData, FilenameMeta})
    fields = map(f -> (f, getfield(x, f)), fieldnames(typeof(x)))
    fields = filter(x -> !isnothing(x[2]), fields)
    map(f -> println(io, f), fields)
end


include("data-utils.jl")
include("data-functions.jl")
include("datamap.jl")
include("helper-functions.jl")
include("diagnostics.jl")

end