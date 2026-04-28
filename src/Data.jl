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
const MODEL_NAME_FIXES = Dict(
    "FGOALS_g2" => "FGOALS-g2", 
    "ACCESS1.3" => "ACCESS1-3",
    "fio-esm" => "FIO-ESM"
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

abstract type AbstractFnFormat end
abstract type FilenameFormat <: AbstractFnFormat end
abstract type ESMVTFormat <: AbstractFnFormat end

struct FF_CMIP <: FilenameFormat end
struct FF_ESMVT_CMIP5 <: FilenameFormat end
struct FF_ESMVT_CMIP6 <: FilenameFormat end
struct FF_ESMVT <: ESMVTFormat end

# Internal conversion
toFF(::Val{:cmip}) = FF_CMIP()
toFF(::Val{:esmvaltool_cmip5})  = FF_ESMVT_CMIP5()
toFF(::Val{:esmvaltool_cmip6}) = FF_ESMVT_CMIP6()
toFF(::Val{:esmvaltool}) = FF_ESMVT()
toFF(f::AbstractFnFormat) = f 
toFF(::Val{f}) where {f} = throw(ArgumentError("invalid filename format :$f, expected one of: :cmip, :esmvaltool, :esmvaltool_cmip5, :esmvaltool_cmip6"))

formatString(::FF_CMIP) = "VARIABLE_TABLEID_MODEL_EXPERIMENT_VARIANT_GRID_TIMERANGE"
formatString(::FF_ESMVT_CMIP5) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE"
formatString(::FF_ESMVT_CMIP6) = "MIP_MODEL_TABLEID_EXPERIMENT_VARIANT_VARIABLE_GRID"

function _parseFormat(format_string::String)
    parts = split(format_string, "_")
    return Dict(Symbol(lowercase(p)) => i for (i, p) in enumerate(parts))
end
# # precomputed at definition time
const CMIP5_FIELD_INDICES = _parseFormat(formatString(FF_ESMVT_CMIP5()))
const CMIP6_FIELD_INDICES = _parseFormat(formatString(FF_ESMVT_CMIP6()))
const CMIP_FIELD_INDICES = _parseFormat(formatString(FF_CMIP()))

struct FilenameMeta
    fn::String
    path::String
    variable::SubString{String}
    tableid::SubString{String}
    model::SubString{String}
    experiment::SubString{String}
    variant::SubString{String}
    grid::SubString{String} # in CMIP standard only given for CMIP6, not CMIP5
    mip::Union{Nothing, SubString{String}} # in CMIP standard only given for CMIP5, not CMIP6
    timerange::Union{Nothing, SubString{String}}# only given in filenames of CMIP standard
    member_id::String
    
    function FilenameMeta(
        fn::String, 
        path::String,
        variable::SubString{String}, 
        tableid::SubString{String}, 
        model::SubString{String}, 
        experiment::SubString{String}, 
        variant::SubString{String}, 
        grid::SubString{String}, 
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

function FilenameMeta(;
    fn::String,
    path::String,
    variable::SubString{String}, 
    tableid::SubString{String}, 
    model::SubString{String}, 
    experiment::SubString{String}, 
    variant::SubString{String}, 
    grid::SubString{String}, 
    mip::Union{Nothing, SubString{String}}, 
    timerange::Union{Nothing, SubString{String}}
)
    FilenameMeta(fn, path, variable, tableid, model, experiment, variant, grid, mip, timerange)
end



const META_FIELDS = fieldnames(FilenameMeta)
const PreviewMap = Dict{String, Vector{FilenameMeta}}

# Constraint has the same fields as FileNameMeta

@kwdef struct Constraint
    fn::Vector{String}
    path::Vector{String}
    variable::Vector{String}
    tableid::Vector{String}
    model::Vector{String}
    experiment::Vector{String}
    variant::Vector{String}
    grid::Vector{String} # in CMIP standard only given for CMIP6, not CMIP5
    mip::Vector{String} # in CMIP standard only given for CMIP5, not CMIP6
    timerange::Vector{String} # only given in filenames of CMIP standard
end


function Constraint(constraint::Dict{Symbol, <:AbstractArray{String}})
    names = fieldnames(FilenameMeta)
    if any(x -> !(x in names), keys(constraint)) 
        throw(ArgumentError("Allowed keys in constraint are: $names"))
    end
    vals = values(constraint)
    if all(isempty, vals)
        throw(ArgumentError("constraint dict is empty!"))
    end
    Constraint(
        fn = get(constraint, :fn, String[]),
        path = get(constraint, :path, String[]),
        variable = get(constraint, :variable, String[]),
        tableid = get(constraint, :tableid, String[]),
        model = get(constraint, :model, String[]),
        experiment = get(constraint, :experiment, String[]),
        variant = get(constraint, :variant, String[]),
        grid = get(constraint, :grid, String[]),
        mip = get(constraint, :mip, String[]),
        timerange = get(constraint, :timerange, String[])
    )
end


# ::MIME"text/plain" : for REPL output
function Base.show(io::IO, ::MIME"text/plain", x::DataMap) #Dict{String, YAXArray})
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

function Base.show(io::IO, ::MIME"text/plain", x::Union{MetaData, FilenameMeta})
    fields = map(f -> (f, getfield(x, f)), fieldnames(typeof(x)))
    fields = filter(field -> !isnothing(field[2]), fields)
    map(f -> println(io, f), fields)
end


include("data-utils.jl")
include("data-functions.jl")
include("helper-functions.jl")
include("diagnostics.jl")



end