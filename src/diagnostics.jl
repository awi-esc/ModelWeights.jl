using DimensionalData
using GLM
using YAXArrays



"""
    computeGlobalMeans(data::AbtractArray)

Return a YAXArray with area-weighted global means for each model in `data`. 

Missing data is accounted for in the area-weights. 

# Arguments:
- `data::YAXArray`: must have dimensions 'lon' and 'lat'.
"""
function computeGlobalMeans(data::YAXArray)
    if !hasdim(data, :lon) || !hasdim(data, :lat)
        msg = "Global means can only be computed for data with :lon, :lat dimensions! Given: $(dims(data))"
        throw(ArgumentError(msg))
    end
    meta = deepcopy(data.properties)
    meta["_statistic"] = hasdim(data, :time) ? "GM-ts" : "GM"
    meta["_id"] = buildMetaDataID(meta)

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
    computeAnomalies(orig_data::YAXArray, ref_data::YAXArray; stats_name::String="ANOM")

Return difference of `orig_data` and `ref_data`.  

The id of the original data and of the reference data is added to the metadata.
"""
function computeAnomalies(orig_data::YAXArray, ref_data::YAXArray; stats_name::String="ANOM")
    dimension, models = getDimsModel(orig_data)
    if !hasdim(ref_data, dimension) || models != Array(dims(ref_data, dimension))
        throw(
            ArgumentError(
                "Original and reference data must contain exactly the same models!",
            ),
        )
    end
    if orig_data.properties["units"] != orig_data.properties["units"]
        @warn "Data and reference data are given in different units! NO ANOMALIES computed!"
        return nothing
    end
    anomalies_metadata = deepcopy(orig_data.properties)
    anomalies_metadata["_statistic"] = stats_name
    anomalies_metadata["_id"] = buildMetaDataID(anomalies_metadata)

    anomalies_metadata["_ref_data_id"] = ref_data.properties["_id"]
    anomalies_metadata["_orig_data_id"] = orig_data.properties["_id"]

    anomalies = @d orig_data .- ref_data
    return YAXArray(dims(orig_data), Array(anomalies), anomalies_metadata)
end

function computeAnomaliesGM(
    orig_data::YAXArray; ref_data::Union{YAXArray, Nothing}=nothing, stats_name::String="ANOM-GM"
)
    gms = isnothing(ref_data) ? computeGlobalMeans(orig_data) : computeGlobalMeans(ref_data)
    return computeAnomalies(orig_data, gms; stats_name)
end


function computeSTD(data::YAXArray, dimension::Symbol)
    meta_new = deepcopy(data.properties)
    meta_new["_statistic"] = "STD"
    meta_new["_id"] = buildMetaDataID(meta_new)
    standard_devs = dropdims(std(data, dims = dimension), dims = dimension)
    return YAXArray(otherdims(data, dimension), standard_devs, meta_new)
end