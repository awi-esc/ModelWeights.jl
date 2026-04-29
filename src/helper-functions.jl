# general helper functions

"""
    average(data::YAXArray, dimension::Symbol)

Wrapper function for mean that removes 'dimension' from the returned YAXArray.
"""
function average(data::YAXArray, dimension::Symbol)
    avg = mean(data; dims=dimension)    
    indices = map(x -> x == dimension ? 1 : Colon(), dimNames(data))
    return avg[indices...]
end

"""
    add(data::YAXArray, dimension::Symbol)

Wrapper function for sum that removes 'dimension' from the returned YAXArray.
"""
function add(data::YAXArray, dimension::Symbol)
    result = sum(data; dims=dimension)    
    indices = map(x -> x == dimension ? 1 : Colon(), dimNames(data))
    return result[indices...]
end


function sharedKeys(d1, d2)
    keys1 = collect(keys(d1))
    keys2 = collect(keys(d2))
    return intersect(keys1, keys2)
end


"""
    throwErrorIfDimMissing(data::YAXArray, dim::Symbol)
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
    if any(ismissing, df)
        base = "Data contains missing values"
        msg = isempty(name) ? base : "$base; $name"
        @warn msg
    end
    return nothing
end


function currentTime(;add_hour::Bool=true)
    current_day = string(Dates.today())
    if add_hour
        current_day = current_day * '_' * Dates.format(Dates.now(), "HH_MM")
    end
    return current_day
end


"""
    individuatePath(target::String)

If file at `target` path exists, concatenate `target` with current timestep, otherwise 
simply return `target`.
"""
function individuatePath(target::String; add_hour::Bool = true)
    target_dir = dirname(target)
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    if isfile(target)
        msg = "$target already exisits, will use "
        target = joinpath(target_dir, join([currentTime(;add_hour), basename(target)], "_"))
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
            model_dim = modelDim(data)
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


function absent(x::AbstractArray)
    isempty(x) || all(absent, x)
end

function absent(x::String)
    isnothing(x) || isempty(x)
end

function absent(x::Dict)
    isempty(x)
end

function absent(x::Nothing) 
    return true
end

"""
    setDim(
        data::YAXArray, dim::Symbol, dim_name::Symbol, dim_vals::AbstractVector{T}
    ) where T

Rename dimension `dim` in `data`to `dim_name` and set dimension values of dimension `dim` to `dim_vals`.

# Arguments:
- `data::YAXArray`:
- `dim::Symbol`: Name of the dimension to be changed.
- `dim_name::Symbol`: New dimension name.
- `dim_vals::AbstractVector{T}`:  New dimension values.
"""
function setDim(
    data::YAXArray, dim::Symbol, dim_name::Symbol, dim_vals::AbstractVector{T}
) where T
    data = setDim(data, dim, dim_name)
    return setDim(data, dim_name, dim_vals)
end


"""
    setDim(data::YAXArray, dim::Symbol, dim_name::Symbol)

Rename dimension `dim` to `dim_name` in `data`.

# Arguments:
- `data::YAXArray`:
- `dim::Symbol`: Name of the dimension to be changed.
- `dim_name::Symbol`: New dimension name.
"""
function setDim(data::YAXArray, dim::Symbol, dim_name::Symbol)
    data = DimensionalData.set(data, dim => dim_name)
    return data
end


"""
    setDim(data::YAXArray, dim::Symbol, dim_vals::AbstractVector{T}) where T

Set dimension values of dimension `dim` in `data` to `dim_vals`.

# Arguments:
- `data::YAXArray`:
- `dim::Symbol`: Name of the dimension to be changed.
- `dim_vals::Vector{String}`:  New dimension values.
"""
function setDim(data::YAXArray, dim::Symbol, dim_vals::AbstractVector{T}) where T
    data = DimensionalData.set(data, dim => dim_vals)
    return data
end


"""
    dimNames(data::YAXArray)

Return the names of the dimensions of `data` as vector of symbols.
"""
function dimNames(data::YAXArray)
    #return collect(DimensionalData.name.(dims(data)))
    return DimensionalData.name.(dims(data))
end


function renameDict!(
    data::Dict{T, V}, ids::AbstractVector{T}, ids_new::AbstractVector{T}
) where {T <: Union{String, Symbol}, V}
    if length(ids) != length(ids_new)
        throw(ArgumentError("The nb of old and new ids must be the same!"))
    end
    if any(id -> !haskey(data, id), ids)
        #@warn "Dictionary does not contain all keys $ids, no key renamed."
        return data
    end
    vals = [data[id] for id in ids]
    for (id, id_new, val) in zip(ids, ids_new, vals)
        delete!(data, id)
        data[id_new] = val
    end
    return data
end


"""
    normalizeDict(data::Dict{String, T}; remove_zero::Bool=true) where T <: Number

Normalize values for every entry in `data` such that they sum up to 1. If remove_zero is
true (default), the returned dictionary does not contain entries for which values were 0.
"""
function normalizeDict(data::Dict{String, T}; remove_zero::Bool=true) where T <: Number
    result = Dict{String, Float64}()
    total = sum(values(data))
    if total == 0
        throw(ArgumentError("Cannot normalize dict by total value of 0!"))
    end
    for (k, v) in data
        if v != 0
            result[k] = data[k] / total
        elseif !remove_zero
            result[k] = 0
        end
    end
    return result
end


"""
    dict2YAX(data::Dict{String, T}) where {T <: Number}

Convert dictionary `data` into a YAXArray with new dimension `dim_name` 
(default is :diagnostic) whose lookup names are the keys of `data`.
"""
function dict2YAX(data::Dict{String, T}; dim_name::Symbol = :diagnostic) where {T <: Number}
    yax = YAXArray(
        (Dim{:diagnostic}(collect(keys(data))),), Array{T}(undef, length(data))
    )
    for k in keys(data)
        yax[diagnostic = At(k)] = data[k]
    end
    if dim_name != :diagnostic
        setDim(yax, :diagnostic, dim_name)
    end
    return yax
end


""" 
    lon360to180(lon::T) where {T <: Real}
    
Convert longitudes measured from 0° to 360° into  -180° to 180° scale.
"""
function lon360to180(lon::T) where {T <: Real}
    return lon > 179 ? lon - 360 : lon
end

function lon360to180(data::YAXArray)
    return setDim(data, :lon, lon360to180.(lookup(data, :lon)))
end


""" 
    lon180to360(lon::T) where {T <: Real}

Convert longitudes measured from -180° to 180° into 0° to 360° scale. For western hemisphere 
(negative longitudes) add 360.
"""
function lon180to360(lon::T) where {T <: Real}
    return ifelse(lon < 0, lon + 360, lon)
end

function lon180to360(data::YAXArray)
    return setDim(data, :lon, lon180to360.(lookup(data, :lon)))
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
function longitudesEastWest(data::AbstractArray{T}) where T <: Number
    data = lon180to360(data)
    longitudes = lookup(data, :lon)
    indices_east = findall(x -> x < 180, longitudes)
    indices_west = findall(x -> x >= 180, longitudes)
    return (east = indices_east, west = indices_west)
end


function countMap(data::AbstractVector{T}) where T    
    counts = Dict{T, Int}()
    for cat in data
        counts[cat] = get(counts, cat, 0) + 1
    end
    return counts
end


function insertSingletonDim(A::AbstractArray, idx::Int)
    old_size = size(A)
    new_size = ntuple(i -> i < idx ? old_size[i] : (i > idx ? old_size[i-1] : 1), ndims(A)+1)
    return reshape(A, new_size)
end

function insertSingletonDim(A::YAXArray, idx::Int, name::Symbol, val::String)
    return _insertSingletonDim(A, idx, Val(name), val)
end

function _insertSingletonDim(A::YAXArray, idx::Int, ::Val{name}, val::String) where {name}
    A_mat = insertSingletonDim(parent(A), idx)
    old_dims = DimensionalData.dims(A)
    new_dim = Dim{name}([val])

    new_dims = DimensionalData.DimTuple((
        old_dims[1:idx-1]...,
        new_dim,
        old_dims[idx:end]...
    ))
    return YAXArray(new_dims, A_mat, copy(A.properties))
end


function dimsEqualButOne(data::AbstractArray{<:YAXArray}; except_dim::Symbol = :time)
    return _dimsEqualButOne(data, Val(except_dim))
end

function _dimsEqualButOne(data::AbstractArray{<:YAXArray}, ::Val{D}) where {D}
    ref = size(otherdims(first(data), D))
    for i in eachindex(data[2:end])
        size(otherdims(data[i], D)) == ref || return false
    end
    return true
end