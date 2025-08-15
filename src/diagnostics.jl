function activeDiagnostics(config::Dict{String, Number})
    return filter(k -> config[k] > 0, collect(keys(config)))
end


"""
    globalMeans(data::YAXArray)

Return a YAXArray with area-weighted global means for each model in `data`. 

Missing data is accounted for in the area-weights. 

# Arguments:
- `data::YAXArray`: must have dimensions 'lon' and 'lat'.
"""
function globalMeans(data::YAXArray)
    throwErrorIfNotLonLat(data)
    latitudes = collect(lookup(data, :lat))
    mask = Bool.(mapslices(model -> ismissing.(model), data; dims=(:lon, :lat)))
    aw_mat = areaWeightMatrix(latitudes, Array(mask)) # is normalized, lon x lat
    # make sure that units are identical across models 
    units = data.properties["units"]
    if length(unique(units)) != 1
        if any(x -> x != "degC" && x != "K", units)
            throw(ArgumentError("Data for all models must be defined in the same units to compute global mean! Only degC and K are converted into degC."))
        end
        data = kelvinToCelsius(data)
    end
    gms = dropdims(mapslices(x -> sum(skipmissing(x)), aw_mat .* Array(data); dims=(1, 2)); dims=(1,2))
    return YAXArray(otherdims(data, (:lon, :lat)), Array(gms), deepcopy(data.properties))
end


"""
    anomalies(orig_data::YAXArray, ref_data::YAXArray)

Compute difference of `orig_data` and `ref_data`.  
"""
function anomalies(orig_data::YAXArray, ref_data::YAXArray)
    data = YAXArray(orig_data.axes, Array(orig_data.data), deepcopy(orig_data.properties))
    dimension = modelDim(data)
    if !hasdim(ref_data, dimension)
        err_msg = "To compute anomalies, ref data must have same model dimension as data (found data: $dimension, found ref_data: $(dims(ref_data)))"
        throw(ArgumentError(err_msg))
    end
    dims_orig = Array(dims(data, dimension))
    dims_ref = Array(dims(ref_data, dimension))
    if length(dims_orig) > length(dims_ref) 
        throw(ArgumentError("To compute anomalies, reference data must contain all models of original data!"))
    end
    if data.properties["units"] != data.properties["units"]
        @warn "Data and reference data are given in different units! NO ANOMALIES computed!"
        return nothing
    end
    if length(dims_orig) < length(dims_ref)
        indices = findall(x -> x in dims_orig, dims_ref)
        ref_data = indexModel(ref_data, (dimension,), indices)
    end
    return @d data .- ref_data
end


function anomaliesGM(orig_data::YAXArray; ref_data::Union{YAXArray, Nothing} = nothing)
    gms = isnothing(ref_data) ? globalMeans(orig_data) : globalMeans(ref_data)
    return anomalies(orig_data, gms)
end


function standardDev(data::YAXArray, dimension::Symbol)
    standard_devs = dropdims(Statistics.std(data, dims = dimension), dims = dimension)
    return YAXArray(otherdims(data, dimension), standard_devs, deepcopy(data.properties))
end