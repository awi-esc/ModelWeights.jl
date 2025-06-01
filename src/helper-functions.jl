using DimensionalData
using YAXArrays


function warn_if_identical_keys(d1, d2)
    keys1 = collect(keys(d1))
    keys2 = collect(keys(d2))
    duplicates = intersect(keys1, keys2)
    if length(duplicates) > 0
        @warn "Merged dictionaries both contain keys: $duplicates. Values from second dictionary are used!"
    end
    return nothing
end


"""
    throwErrorIfModelDimMissing(data::AbstractArray)
"""
function throwErrorIfModelDimMissing(data::AbstractArray)
    if !hasdim(data, :model) && !hasdim(data, :member)
        throw(ArgumentError("Data must have dimension 'model' or 'member'."))
    end
    return nothing
end


"""
    getDimsModel(data::AbstractArray)

Return Tuple with the model dimension of `data` ('model' or 'member')
as Symbol at position 1 and the respetive dimension names at position 2.
"""
function getDimsModel(data::AbstractArray)
    throwErrorIfModelDimMissing(data)
    dim_symbol = hasdim(data, :model) ? :model : :member
    return (dim_symbol, Array(dims(data, dim_symbol)))
end


"""
    getAtModel(data::AbstractArray, dimension::Symbol, model::String)

Return `data` where `dimension` (member or model) has value `model`.
"""
function getAtModel(data::AbstractArray, dimension::Symbol, model::String)
    throwErrorIfModelDimMissing(data)
    return dimension == :model ? data[model=At(model)] : data[member=At(model)]
end


"""
    getByIdxModel(data::AbstractArray, dimension::Symbol, indices::Vector)

Return `data` where `dimension` (member or model) has value `model`.
"""
function getByIdxModel(data::AbstractArray, dimension::Symbol, indices::Vector)
    throwErrorIfModelDimMissing(data)
    return dimension == :model ? data[model=indices] : data[member=indices]
end


""" 
    putAtModel!(data::AbstractArray, dimension::Symbol, model::String, input)
"""
function putAtModel!(data::AbstractArray, dimension::Symbol, model::String, input)
    throwErrorIfModelDimMissing(data)
    if dimension == :model
        data[model=At(model)] = input
    else
        data[member=At(model)] = input
    end
    return nothing
end


"""
    mergeMetaDataPaths(meta1::Dict{String, Any}, meta2::Dict{String, Any})

Return vector of unique paths of mergee set of paths from `meta1` and `meta2`.
"""
function mergeMetaDataPaths(meta1::Dict{String,Any}, meta2::Dict{String,Any})
    paths = copy(meta1["_paths"])
    append!(paths, meta2["_paths"])
    return unique(paths)
end


function buildMetaDataID(variable::String, statistic::String, alias::String)
    return join([variable, statistic, alias], "_")
end

function buildMetaDataID(meta::Dict{String,Any})
    return join([meta["_variable"], meta["_statistic"], meta["_alias"]], "_")
end


"""
    getAllCombinations(v...)

Generate vector of strings with all possible combinations of input vectors,
where each combination consists of one element from each input vector, 
concatenated as a string with underscores separating the elements.

# Arguments
- `v...`: A variable number of input vectors.

# Example
```jldoctest
julia> ModelWeights.getAllCombinations(["tos", "tas"], ["CLIM"])
2-element Vector{String}:
 "tos_CLIM"
 "tas_CLIM"
```
"""
function getAllCombinations(v...)
    if any(isempty.(v))
        @warn "At least one input vector is empty -> empty vector returned!"
    end
    combis = Vector{String}()
    for elems in Iterators.product(v...)
        push!(combis, join(elems, "_"))
    end
    return combis
end

function warnIfhasMissing(df::YAXArray; name::String="")
    base = "Data contains missing values"
    msg = isempty(name) ? base : base * "; $(name)"
    if any(ismissing.(df))
        @warn msg
    end
    return nothing
end