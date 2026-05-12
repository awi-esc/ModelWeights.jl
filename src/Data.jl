module Data

export Level, MODEL_MEMBER_DELIM
export DataMap

export subsetModelData, sharedModels
export summarizeMembers, summarizeMembers!
export alignPhysics, addMasks!
export defineDataMap, loadPreprocData

import StatsBase
import Missings:allowmissing
import LinearRegression as lr

using CSV
using Dates
using DataFrames
using DimensionalData
using Distributed
using NCDatasets
using NetCDF
using JLD2
using Serialization
using Setfield
using Statistics
using YAML
using YAXArrays

const MODEL_MEMBER_DELIM = "#"
# maps from file name as retrieved from metadata to name in the filename itself
const MODEL_NAME_FIXES = Dict(
    "FGOALS_g2" => "FGOALS-g2", 
    "ACCESS1.3" => "ACCESS1-3",
    "FIO-ESM" => "fio-esm"
    #"fio-esm" => "FIO-ESM"
)

abstract type AbstractLevel end
abstract type Level <: AbstractLevel end
abstract type NoLevel <: AbstractLevel end

struct LevMember <: Level end
struct LevModel <: Level end
struct LevNone <: NoLevel end

# Internal conversion
toLevel(::Val{:member}) = LevMember()
toLevel(::Val{:model})  = LevModel()
toLevel(::Val{:none}) = LevNone()
toLevel(l::Level) = l 
toLevel(::Val{l}) where {l} = throw(ArgumentError("invalid level :$l, expected one of: :member, :model, :none"))


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

function defineDataMap(data::Vector{<:YAXArray}, ids::Vector{String}; sorted::Bool = true)
    if length(data) != length(ids)
        throw(ArgumentError("data and ids must have same size to build up DataMap"))
    end
    dm = DataMap()
    for (d, id) in zip(data, ids)
        if haskey(dm, id) # merge on model dimension
            model_dim = modelDim(d)
            if modelDim(dm[id]) != model_dim
                throw(ArgumentError("Cant join two identical ids, with :member and :model dimension."))
            end
            dm[id] = mergeYAX(dm[id], d, model_dim)
        else
            dm[id] = d
        end
    end
    return dm
end

abstract type AbstractFnFormat end
abstract type FilenameFormat <: AbstractFnFormat end
abstract type ESMVTFormat <: AbstractFnFormat end

struct FF_CMIP <: FilenameFormat end

struct FF_ESMVT_CMIP5 <: FilenameFormat end
struct FF_ESMVT_CMIP5_TR <: FilenameFormat end
struct FF_ESMVT_CMIP6 <: FilenameFormat end
struct FF_ESMVT_CMIP6_TR <: FilenameFormat end

struct FF_ESMVT_OBS <: FilenameFormat end
struct FF_ESMVT <: ESMVTFormat end

# Internal conversion
toFF(::Val{:cmip}) = FF_CMIP()
toFF(::Val{:esmvaltool_cmip5})  = FF_ESMVT_CMIP5()
toFF(::Val{:esmvaltool_cmip6}) = FF_ESMVT_CMIP6()
toFF(::Val{:esmvaltool}) = FF_ESMVT()
toFF(::Val{:esmvaltool_obs}) = FF_ESMVT_OBS()
toFF(f::AbstractFnFormat) = f 
toFF(::Val{f}) where {f} = throw(ArgumentError("invalid filename format :$f, expected one of: :cmip, :esmvaltool, :esmvaltool_cmip5, :esmvaltool_cmip6, :esmvaltool_obs"))

formatString(::FF_CMIP) = "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT_GRID_TIMERANGE"

formatString(::FF_ESMVT_CMIP5) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE"
formatString(::FF_ESMVT_CMIP5_TR) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_TIMERANGE"

formatString(::FF_ESMVT_CMIP6) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_GRID"
formatString(::FF_ESMVT_CMIP6_TR) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_GRID_TIMERANGE"

formatString(::FF_ESMVT_OBS) = "GRID_MODEL_TYPE_VERSION_TABLEID_VARIABLE"

function _parseFormat(format_string::String)
    parts = split(format_string, "_")
    return Dict(Symbol(lowercase(p)) => i for (i, p) in enumerate(parts))
end
# # precomputed at definition time
const CMIP5_FIELD_INDICES = _parseFormat(formatString(FF_ESMVT_CMIP5()))
const CMIP5_FIELD_INDICES_TR = _parseFormat(formatString(FF_ESMVT_CMIP5_TR()))

const CMIP6_FIELD_INDICES = _parseFormat(formatString(FF_ESMVT_CMIP6()))
const CMIP6_FIELD_INDICES_TR = _parseFormat(formatString(FF_ESMVT_CMIP6_TR()))

const CMIP_FIELD_INDICES = _parseFormat(formatString(FF_CMIP()))
const OBS_FIELD_INDICES = _parseFormat(formatString(FF_ESMVT_OBS()))


abstract type AbstractMeta end
struct ModelMeta <: AbstractMeta
    fn::String
    path::String
    variable::SubString{String}
    tableid::SubString{String}
    model::SubString{String}
    experiment::SubString{String}
    variant::SubString{String}
    grid::Union{Nothing, SubString{String}} # in CMIP standard only given for CMIP6, not CMIP5
    mip::Union{Nothing, SubString{String}} # in CMIP standard only given for CMIP5, not CMIP6
    timerange::Union{Nothing, SubString{String}}# only given in filenames of CMIP standard
    member_id::String
    
    function ModelMeta(
        fn::String, 
        path::String,
        variable::SubString{String}, 
        tableid::SubString{String}, 
        model::SubString{String}, 
        experiment::SubString{String}, 
        variant::SubString{String}, 
        grid::Union{Nothing, SubString{String}}, 
        mip::Union{Nothing, SubString{String}}, 
        timerange::Union{Nothing, SubString{String}}
    )
        member_id = string(model, MODEL_MEMBER_DELIM, variant)
        if !isnothing(mip) && lowercase(mip) == "cmip6"
            member_id = string(member_id, "_", grid)
        end
        new(fn, path, variable, tableid, model, experiment, variant, grid, mip, timerange, member_id)
    end
end

function ModelMeta(;
    fn::String,
    path::String,
    variable::SubString{String}, 
    tableid::SubString{String}, 
    model::SubString{String}, 
    experiment::SubString{String}, 
    variant::SubString{String}, 
    grid::Union{Nothing, SubString{String}}, 
    mip::Union{Nothing, SubString{String}}, 
    timerange::Union{Nothing, SubString{String}}
)
    ModelMeta(fn, path, variable, tableid, model, experiment, variant, grid, mip, timerange)
end

struct ObsMeta <: AbstractMeta
    fn::String
    path::String
    variable::SubString{String}
    tableid::SubString{String} 
    model::SubString{String}
    grid::SubString{String}
    type::SubString{String} 
    version::SubString{String} 
end

function ObsMeta(;
    fn::String, 
    path::String,
    variable::SubString{String},
    tableid::SubString{String}, 
    model::SubString{String}, 
    grid::SubString{String}, 
    type::SubString{String}, 
    version::SubString{String}
)
    ObsMeta(fn, path, variable, tableid, model, grid, type, version)
end

const PreviewMap = Dict{String, Vector{AbstractMeta}}

@kwdef struct Constraint
    filenames::Vector{String}
    variables::Vector{String}
    tableids::Vector{String}
    models::Vector{String}
    experiments::Vector{String}
    variants::Vector{String}
    grids::Vector{String} # in CMIP standard only given for CMIP6, not CMIP5
    mips::Vector{String} # in CMIP standard only given for CMIP5, not CMIP6
    members::Vector{String} # only applicable to model data, not observations
end

const CONSTRAINT_FIELDS = fieldnames(Constraint)

function Constraint(constraint::Dict{Symbol, Vector{String}}; warn::Bool=true)
    all_keys = keys(constraint)
    irrelevant_keys = filter(x -> !(x in CONSTRAINT_FIELDS), all_keys)
    if !isempty(irrelevant_keys) && warn
        @warn "The following entries in constraints are not valid and are not taken into account: $irrelevant_keys.\nValid keys are: $CONSTRAINT_FIELDS"
    end
    Constraint(
        filenames = get(constraint, :filenames, String[]),
        variables = get(constraint, :variables, String[]),
        tableids = get(constraint, :tableids, String[]),
        models = get(constraint, :models, String[]),
        experiments = get(constraint, :experiments, String[]),
        variants = get(constraint, :variants, String[]),
        grids = get(constraint, :grids, String[]),
        mips = get(constraint, :mips, String[]),
        members = get(constraint, :members, String[])
    )
end


const CONSTRAINTS_TO_META = Dict(
    :filenames => :fn,
    :variables => :variable,
    :tableids => :tableid,
    :models => :model,
    :experiments => :experiment,
    :variants => :variant,
    :grids => :grid, 
    :mips => :mip
)



# ::MIME"text/plain" : for REPL output
function Base.show(io::IO, ::MIME"text/plain", x::DataMap)
    println(io, "::DataMap")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::PreviewMap)
    println(io, "::PreviewMap")
    for (k, v) in x
        println(io, "$k: $(size(v))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::Union{MetaData, AbstractMeta, Constraint})
    fields = map(f -> (f, getfield(x, f)), fieldnames(typeof(x)))
    fields = filter(field -> !isnothing(field[2]), fields)
    map(f -> println(io, f), fields)
end


include("data-utils.jl")
include("data-functions.jl")
include("helper-functions.jl")
include("diagnostics.jl")



end