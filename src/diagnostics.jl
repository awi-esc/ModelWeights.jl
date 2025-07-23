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
    throwErrorIfDimMissing(data, [:lon, :lat]; include = :any)
    longitudes = Array(dims(data, :lon))
    latitudes = Array(dims(data, :lat))
    area_weights = approxAreaWeights(coalesce.(latitudes, NaN))
    s = otherdims(data, (:lon, :lat))
    area_weighted_mat =
        isempty(s) ? repeat(area_weights', length(longitudes), 1) :
        repeat(area_weights', length(longitudes), 1, size(s)...)
    masks = ismissing.(data)
    area_weighted_mat = ifelse.(masks .== 1, missing, area_weighted_mat)

    weighted_unnormalized_vals = area_weighted_mat .* data
    normalization = mapslices(area_weighted_mat, dims = ("lon", "lat")) do x
        sum(skipmissing(x))
    end
    weighted_normalized_vals = @d weighted_unnormalized_vals ./ normalization
    gms = mapslices(weighted_normalized_vals, dims = ("lon", "lat")) do x
        sum(skipmissing(x))
    end
    return YAXArray(otherdims(data, (:lon, :lat)), Array(gms), deepcopy(data.properties))
end


"""
    anomalies(orig_data::YAXArray, ref_data::YAXArray)

Compute difference of `orig_data` and `ref_data`.  
"""
function anomalies(orig_data::YAXArray, ref_data::YAXArray)
    dimension = modelDim(orig_data)
    if !hasdim(ref_data, dimension)
        err_msg = "To compute anomalies, ref data must have same model dimension as data (found data: $dimension, found ref_data: $(dims(ref_data)))"
        throw(ArgumentError(err_msg))
    end
    dims_orig = Array(dims(orig_data, dimension))
    dims_ref = Array(dims(ref_data, dimension))
    if length(dims_orig) > length(dims_ref) 
        throw(ArgumentError("To compute anomalies, reference data must contain all models of original data!"))
    end
    if orig_data.properties["units"] != orig_data.properties["units"]
        @warn "Data and reference data are given in different units! NO ANOMALIES computed!"
        return nothing
    end
    if length(dims_orig) < length(dims_ref)
        indices = findall(x -> x in dims_orig, dims_ref)
        ref_data = indexModel(ref_data, (dimension,), indices)
    end
    anomalies = @d orig_data .- ref_data
    return YAXArray(dims(orig_data), anomalies.data, deepcopy(orig_data.properties))
end


function anomaliesGM(orig_data::YAXArray; ref_data::Union{YAXArray, Nothing} = nothing)
    gms = isnothing(ref_data) ? globalMeans(orig_data) : globalMeans(ref_data)
    return anomalies(orig_data, gms)
end


function standardDev(data::YAXArray, dimension::Symbol)
    standard_devs = dropdims(Statistics.std(data, dims = dimension), dims = dimension)
    return YAXArray(otherdims(data, dimension), standard_devs, deepcopy(data.properties))
end