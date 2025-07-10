"""
    globalMeans(data::AbtractArray)

Return a YAXArray with area-weighted global means for each model in `data`. 

Missing data is accounted for in the area-weights. 

# Arguments:
- `data::YAXArray`: must have dimensions 'lon' and 'lat'.
"""
function globalMeans(data::YAXArray)
    if !hasdim(data, :lon) || !hasdim(data, :lat)
        msg = "Global means can only be computed for data with :lon, :lat dimensions! Given: $(dims(data))"
        throw(ArgumentError(msg))
    end
    meta = deepcopy(data.properties)
    meta["statistic"] = hasdim(data, :time) ? "GM-ts" : "GM"
    meta["id"] = buildMetaDataID(meta)

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
    return YAXArray(otherdims(data, (:lon, :lat)), Array(gms), meta)
end


"""
    anomalies(orig_data::YAXArray, ref_data::YAXArray; stats_name::String="ANOM")

Compute difference of `orig_data` and `ref_data`.  

The id of the original data and of the reference data is added to the metadata.
"""
function anomalies(orig_data::YAXArray, ref_data::YAXArray; stats_name::String="ANOM")
    dimension = modelDim(orig_data)
    if !hasdim(ref_data, dimension) || (Array(dims(orig_data, dimension))) != Array(dims(ref_data, dimension))
        err_msg = "To compute anomalies, original and reference data must contain exactly the same models!"
        throw(ArgumentError(err_msg))
    end
    if orig_data.properties["units"] != orig_data.properties["units"]
        @warn "Data and reference data are given in different units! NO ANOMALIES computed!"
        return nothing
    end
    anomalies_metadata = deepcopy(orig_data.properties)
    anomalies_metadata["statistic"] = stats_name
    anomalies_metadata["id"] = buildMetaDataID(anomalies_metadata)

    anomalies_metadata["_ref_data_id"] = ref_data.properties["id"]
    anomalies_metadata["_orig_data_id"] = orig_data.properties["id"]

    anomalies = @d orig_data .- ref_data
    return YAXArray(dims(orig_data), Array(anomalies), anomalies_metadata)
end


function anomaliesGM(
    orig_data::YAXArray; ref_data::Union{YAXArray, Nothing}=nothing, stats_name::String="ANOM-GM"
)
    gms = isnothing(ref_data) ? globalMeans(orig_data) : globalMeans(ref_data)
    return anomalies(orig_data, gms; stats_name)
end


function standardDev(data::YAXArray, dimension::Symbol)
    meta_new = deepcopy(data.properties)
    meta_new["statistic"] = "STD"
    meta_new["id"] = buildMetaDataID(meta_new)
    standard_devs = dropdims(Statistics.std(data, dims = dimension), dims = dimension)
    return YAXArray(otherdims(data, dimension), standard_devs, meta_new)
end