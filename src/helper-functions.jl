# general helper functions

function sharedKeys(d1, d2)
    keys1 = collect(keys(d1))
    keys2 = collect(keys(d2))
    return intersect(keys1, keys2)
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
    mergeMetaDataPaths(meta1::MetaData, meta2::MetaData)

Return vector of unique paths of merged set of paths from `meta1` and `meta2`.
"""
function mergeMetaDataPaths(meta1::MetaData, meta2::MetaData)
    paths = copy(meta1.paths)
    append!(paths, meta2.paths)
    return unique(paths)
end


function buildMetaDataID(meta::MetaData)
    subdir = meta.subdir
    if isempty(subdir)
        subdir = (!isempty(meta.statistic) && !isempty(meta.variable)) ? 
            join([meta.variable, meta.statistic], "_") :
            throw(ArgumentError("Neither subdir nor statistic and variable are in metadata!"))
    end
    if isempty(meta.alias)
        throw(ArgumentError("alias (name of diagnostic in ESMValTool) cannot be empty in MetaData!"))
    end
    return join([subdir, meta.alias], "_")
end

function buildMetaDataID(meta::Dict)
    fn_err(k::String) = throw(ArgumentError("$k not defined in metadata!"))
    subdir = get(meta, "subdir", nothing)
    if isnothing(subdir)
        subdir = join(
            [get(() -> fn_err("subdir, variable"), meta, "variable"),
             get(() -> fn_err("subdir, statistic"), meta, "statistic")], 
             "_"
        )
    end
    return join([subdir, get(() -> fn_err("alias"), meta, "alias")], "_")
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
    return reduce(vcat, map(elems -> join(elems, sep), Iterators.product(v...)))
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

