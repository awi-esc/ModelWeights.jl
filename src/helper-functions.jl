# general helper functions

function sharedKeys(d1, d2)
    keys1 = collect(keys(d1))
    keys2 = collect(keys(d2))
    return intersect(keys1, keys2)
end


"""
    throwErrorIfDimMissing(data::YAXArray)
"""
function throwErrorIfDimMissing(data::YAXArray, dim::Symbol)
    if !hasdim(data, dim)
        throw(ArgumentError("Data must have dimension $(dim). Found: $(name.(dims(data)))"))
    end
    return nothing
end

"""
    throwErrorIfDimMissing(data::YAXArray, dimensions::Vector{Symbol}; include::Symbol=:all)

Throw ArgumentError if `data` does not have ALL dimensions in `dimensions` when 
`include=:all` (default), or if `data` does not have ANY dimension in `dimensions`.
"""
function throwErrorIfDimMissing(data::YAXArray, dimensions::Vector{Symbol}; include::Symbol=:all)
    if include == :all && any(x -> !hasdim(data, x), dimensions)
        throw(ArgumentError("Data must have all of dimensions $(dimensions). Found: $(DimensionalData.name.(DimensionalData.dims(data)))"))
    elseif  all(x -> !hasdim(data, x), dimensions)
        throw(ArgumentError("Data must have any of dimensions $(dimensions). Found: $(DimensionalData.name.(DimensionalData.dims(data)))"))
    end
    return nothing
end


function throwErrorIfNotLonLat(data::YAXArray)
    try
        dimensions = DimensionalData.name.(dims(data))
        if dimensions[1:2]  != (:lon, :lat)
            throw(ArgumentError("Input data must have :lon as first and :lat as second dimension! Found:$(dimensions)"))
        end
    catch _
        throw(ArgumentError("Input data must have :lon as first and :lat as second dimension! Found:$(DimensionalData.name.(dims(data)))"))
    end
    return nothing
end

"""
    combineAll(v::Vararg{Vector{String}}; sep::String="_")

Generate vector with all possible combinations of strings, each of which is a concatenation 
of elements from each input vector in the given order, concatenated by `sep` (default: _).

# Arguments
- `v::Vararg{Vector{String}}`: A variable number of input vectors.
- `sep::String`: seperator

# Example
```jldoctest
julia> ModelWeights.Data.combineAll(["tos", "tas"], ["CLIM"])
2-element Vector{String}:
 "tos_CLIM"
 "tas_CLIM"
```
"""
function combineAll(v::Vararg{Vector{String}}; sep::String="_")
    if any(isempty.(v))
        @warn "At least one input vector is empty -> empty vector returned!"
        return Vector{String}()
    end
    combinations =  map(elems -> join(elems, sep), Iterators.product(v...))
    # init in reduce is important! otherwise if combinations is just one element, that 
    # element would be returned as String!
    return reduce(vcat, combinations; init=Vector{String}())
end


function warnIfhasMissing(df::YAXArray; name::String="")
    base = "Data contains missing values"
    msg = isempty(name) ? base : base * "; $(name)"
    if any(ismissing.(df))
        @warn msg
    end
    return nothing
end


function currentTime()
    currentDay = string(Dates.today()) * '_'
    return currentDay * Dates.format(Dates.now(), "HH_MM")
end


"""
    individuatePath(target::String)

If file at `target` path exists, concatenate `target` with current timestep, otherwise 
simply return `target`.
"""
function individuatePath(target::String)
    target_dir = dirname(target)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    if isfile(target)
        msg = "$target already exisits, will use "
        target = joinpath(target_dir, join([currentTime(), basename(target)], "_"))
        @info msg * "$target instead."
    end
    return target
end


"""
    kelvinToCelsius(data::YAXArray)

Return a copy of `data` with values given in Kelvin covnerted into Degree Celsius.

"""
function kelvinToCelsius(data::YAXArray)
    df = YAXArray(data.axes, copy(data.data), deepcopy(data.properties))
    kelvinToCelsius!(df)
    return df
end

function kelvinToCelsius!(data::YAXArray)
    units = data.properties["units"]
    if isa(units, String) && units == "K"
        data = data .- 273.15
        data.properties["units"] = "degC"
    elseif isa(units, Vector)
        indices = findall(units .== "K")
        if !isempty(indices)
            model_dim = modelDims(data)
            if model_dim == :member
                data[member = indices] .= data[member = indices] .- 273.15
            else
                data[model = indices] .= data[model = indices] .- 273.15
            end
            units[indices] .= "degC"
        end
    end
    return nothing
end



function absent(x::Union{Vector, String, Dict, Nothing}) 
    return isnothing(x) || isempty(x)
end

"""
    setDim(
        data::YAXArray, 
        dim::Union{String, Symbol},
        dim_name::Union{Nothing, String, Symbol},
        dim_vals::Union{Nothing, Vector{String}}
    )

Rename dimension `dim` to `dim_name` and/or its values to `dim_vals`.
If `dim_name` or `dim_vals` is nothing, name/vals don't change.

# Arguments:
- `data::YAXArray`:
- `dim::Union{String, Symbol}`: Name of the dimension to be changed.
- `dim_name::Union{Nothing, String, Symbol}`: New dimension name.
- `dim_vals::Union{Nothing, Vector{String}}`:  New dimension values.
"""
function setDim(
    data::YAXArray, 
    dim::Union{String, Symbol},
    dim_name::Union{Nothing, String, Symbol},
    dim_vals::Union{Nothing, Vector{T}}
) where T <: Union{String, Number}
    if !isnothing(dim_vals)
        data = DimensionalData.set(data, Symbol(dim) => dim_vals)
    end
    if !isnothing(dim_name)
        data = DimensionalData.set(data, Symbol(dim) => Symbol(dim_name))
    end
    return data
end

"""
    dimNames(data::YAXArray)

Return the names of the dimensions of `data` as vector of symbols.
"""
function dimNames(data::YAXArray)
    return collect(DimensionalData.name.(dims(data)))
end


function renameDict!(
    data::Dict{T, V}, ids::AbstractVector{T}, ids_new::AbstractVector{T}
) where {T <: Union{String, Symbol}, V}
    if any(id -> !(id in keys(data)), ids)
        @warn "Dictionary does not contain all keys $ids, no key renamed."
        return nothing
    end
    if length(ids) != length(ids_new)
        throw(ArgumentError("The nb of old and new ids must be the same!"))
    end
    map(ids, ids_new) do id, id_new
        data[id_new] = data[id]
        delete!(data, id)
    end
    return nothing
end


"""
    normalizeDict(data::Dict{String, <:Number})

Normalize values for every entry in `data` such that they sum up to 1. If remove_zero is
true (default), the returned dictionary does not contain entries for which values were 0.
"""
function normalizeDict(data::Dict{String, <:Number}; remove_zero::Bool=true)
    result = Dict{String, Float64}()
    total = sum(values(data))
    data = remove_zero ? filter(((k, v),) -> v != 0, data) : data
    for k in keys(data)
        result[k] = data[k] ./ total 
    end
    return result
end


"""
    dict2YAX(data::Dict{String, <:Number})

Convert dictionary `data` into a YAXArray with new dimension `dim_name` 
(default is :diagnostic) whose lookup names are the keys of `data`.
"""
function dict2YAX(data::Dict{String, <:Number}; dim_name::Symbol = :diagnostic)
    yax = YAXArray(
        (Dim{:diagnostic}(collect(keys(data))),), Array{Float64}(undef, length(data))
    )
    for k in keys(data)
        yax[diagnostic = At(k)] = data[k]
    end
    if dim_name != :diagnostic
        setDim(yax, :diagnostic, dim_name, nothing)
    end
    return yax
end


""" 
    lon360to180(lon::Number)
    
Convert longitudes measured from 0° to 360° into  -180° to 180° scale.
"""
function lon360to180(lon::Number)
    return lon > 179 ? lon - 360 : lon
end

function lon360to180(data::YAXArray)
    return setDim(data, :lon, nothing, lon360to180.(lookup(data, :lon)))
end


""" 
    lon180to360(lon::Number)

Convert longitudes measured from -180° to 180° into 0° to 360° scale. For western hemisphere 
(negative longitudes) add 360.
"""
function lon180to360(lon::Number)
    return ifelse(lon < 0, lon + 360, lon)
end

function lon180to360(data::YAXArray)
    return setDim(data, :lon, nothing, lon180to360.(lookup(data, :lon)))
end


"""
    sortLongitudesWest2East(data::AbstractArray)

Arrange 'data' such that western latitudes come first, then eastern latitudes.
"""
function sortLongitudesEast2West(data::AbstractArray)
    indices = longitudesEastWest(data)
    return data[lon = vcat(indices.east, indices.west)]
end

"""
    sortLongitudesWest2East(data::AbstractArray)

Arrange 'data' such that western latitudes come first, then eastern latitudes.
"""
function sortLongitudesWest2East(data::AbstractArray)
    indices = longitudesEastWest(data)
    # east = longitudes[longitudes .< 180]
    # west = longitudes[longitudes .>= 180]
    # sorted_lon = vcat(west, east)

    # # necessary to specify that lookup dimension values aren't sorted anymore!
    # # otherwise Selector At won't work! does seem to work don't know TODO: CHECK THIS!!
    # lookup_lon = Lookups.Sampled(
    #     sorted_lon;
    #     span = Lookups.Irregular(minimum(longitudes), maximum(longitudes)),
    #     order = Lookups.ForwardOrdered(),
    # )
    # data = data[lon = At(sorted_lon)]
    # data = DimensionalData.Lookups.set(data, lon = lookup_lon)
    return data[lon = vcat(indices.west, indices.east)]
end


"""
    sortLongitudesWest2East(data::AbstractArray)

Arrange 'data' such that western latitudes come first, then eastern latitudes.

Return NamedTuple with fields 'east' and 'west' pointing to vectors of indices for sorted longitudes.
"""
function longitudesEastWest(data::AbstractArray)
    data = lon180to360(data)
    longitudes = lookup(data, :lon)
    indices_east = findall(x -> x < 180, longitudes)
    indices_west = findall(x -> x >= 180, longitudes)
    return (east = indices_east, west = indices_west)
end


function countMap(data::Vector{T}) where T <: Any
    counts = Dict()
    for cat in unique(data)
        nb = count(i -> i == cat, data)
        counts[cat] = nb
    end
    return counts
end