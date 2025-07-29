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
        throw(ArgumentError("Data must have dimension $(dim). Found: $(dims(data))"))
    end
    return nothing
end

"""
    throwErrorIfDimMissing(data::YAXArray, dims::Vector{Symbol}; include::Symbol=:all)

Throw ArgumentError if `data` does not have ALL dimensions in `dims` when `include=:all` 
(default), or if `data` does not have ANY dimension in `dims`.
"""
function throwErrorIfDimMissing(data::YAXArray, dims::Vector{Symbol}; include::Symbol=:all)
    if include == :all && any(x -> !hasdim(data, x), dims)
        throw(ArgumentError("Data must have all of dimensions $(dims). Found: $(dims(data))"))
    elseif  all(x -> !hasdim(data, x), dims)
        throw(ArgumentError("Data must have any of dimensions $(dims). Found: $(dims(data))"))
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
    kelvinToCelsius(data::AbstractArray)

Return a copy of `data` with values given in Kelvin covnerted into Degree Celsius.

"""
function kelvinToCelsius(data::YAXArray)
    units = data.properties["units"]
    df = deepcopy(data)
    if isa(units, String) && units == "K"
        df = df .- 273.15
        df.properties["units"] = "degC"
    elseif isa(units, Vector)
        indices = findall(units .== "K")
        if !isempty(indices)
            if hasdim(df, :member)
                df[member=indices] .= df[member=indices] .- 273.15
            else
                df[model=indices] .= df[model=indices] .- 273.15
            end
            df.properties["units"] = "degC"
        end
    end
    return df
end


"""
    kelvinToCelsius!(datamap::DataMap)

Modify entries of `datamap` s.t. all data is given in Degree Celsius (instead) of Kelvin.
"""
function kelvinToCelsius!(datamap::DataMap)
    for (id, da) in datamap
        datamap[id] = kelvinToCelsius(da)
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
    dim_vals::Union{Nothing, Vector{String}}
)
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
    return map(d -> typeof(d).parameters[1], dims(data))
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