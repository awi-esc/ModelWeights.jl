using DimensionalData
using GLM
using YAXArrays


"""
    addDiagnostic!(datamap::DataMap, fn::Function, args...; kwargs...)

Compute a diagnostic by running `fn` with positional arguments `args` and 
keyword arguments `kwargs` and add result to `datamap` at the id of the 
computed result.
"""
function addDiagnostic!(
    datamap::DataMap, fn::Function, id::String, args...; kwargs...
)
    @info "run $(String(Symbol(fn))) for $id..."
    data = fn(datamap[id], args...; kwargs...)
    datamap[data.properties["_id"]] = data
    return nothing
end


"""
    computeLinearTrend(data::YAXArray; full_predictions::Bool=true)

Compute linear trend as ordinary least squares for timeseries data.

# Arguments:
- `data::YAXArray`: must have dimension :time.
- `full_predictions::Bool`: if false, return only slope, otherwise return 
predicted value.
"""
function computeLinearTrend(data::YAXArray; full_predictions::Bool=false)
    if !hasdim(data, :time)
        msg = "Linear trend can only be computed for timeseries data (with dimension :time)!"
        throw(ArgumentError(msg))
    end
    x = Dates.year.(Array(data.time))
    meta = deepcopy(data.properties)
    stats = full_predictions ? "TREND-pred" : "TREND"
    meta["_id"] = replace(meta["_id"], meta["_statistic"] => stats)
    meta["_statistic"] = stats
 
    function fn(y)
        y = Array(y)
        indices = findall(.!(ismissing.(y)))
        # TODO: this is slow and should be changed!
        ols = lm(@formula(Y~X), DataFrame(X=x[indices], Y=Float64.(y[indices])))
        if full_predictions
            y_hat = Array{Union{Missing, Float64}}(undef, size(y))
            y_hat[indices] =  predict(ols)
        else
            y_hat = coef(ols)[2]
        end
        return y_hat
    end
    trends = mapslices(fn, data; dims = (:time,))
    return YAXArray(dims(trends), Array(trends), meta) 
end


"""
    addLinearTrend!(
        data::DataMap;     
        statistic::String="CLIM-ann", 
        ids_ts::Vector{String}=Vector{String}(), 
        full_predictions::Bool=true
    )

Add computed linear trend (computeLinearTrend) for all datasets in `data` with an 
id in `ids_ts`, or if `ids_ts` is not given for all datasets with `statistic` 
(as specified in metadata _statistic), where the default is the annual 
climatologies (CLIM-ann).
"""
function addLinearTrend!(
    data::DataMap; 
    statistic::String="CLIM-ann", 
    ids_ts::Vector{String}=Vector{String}(), 
    full_predictions::Bool=true
)
    if isempty(ids_ts)
        ids_ts = filter(
            id -> data[id].properties["_statistic"] == statistic, keys(data)
        ) 
    end
    for id in ids_ts
        addDiagnostic!(data, computeLinearTrend, id; full_predictions)
    end
    return nothing
end


"""
    computeGlobalMeans(data::AbtractArray)

Return a YAXArray with area-weighted global means for each model in `data`. 

Missing data is accounted for in the area-weights. 

# Arguments:
- `data::YAXArray`: must have dimensions 'lon' and 'lat'.
"""
function computeGlobalMeans(data::YAXArray)
    longitudes = Array(dims(data, :lon))
    latitudes = Array(dims(data, :lat))
    meta = makeMetadataGMS(data.properties, hasdim(data, :time))
    
    area_weights = approxAreaWeights(coalesce.(latitudes, NaN))
    s = otherdims(data, (:lon, :lat))
    area_weighted_mat = isempty(s) ? 
        repeat(area_weights', length(longitudes), 1) :
        repeat(area_weights', length(longitudes), 1, size(s)...)
    masks = ismissing.(data)
    area_weighted_mat = ifelse.(masks .== 1, missing, area_weighted_mat)

    weighted_unnormalized_vals = area_weighted_mat .* data
    normalization = mapslices(area_weighted_mat, dims=("lon", "lat")) do x 
        sum(skipmissing(x))
    end
    weighted_normalized_vals = @d weighted_unnormalized_vals ./ normalization
    gms = mapslices(weighted_normalized_vals, dims=("lon", "lat")) do x
        sum(skipmissing(x))
    end 
    return YAXArray(otherdims(data, (:lon, :lat)), Array(gms), meta)
end


"""
    addGlobalMeans!(data::DataMap; ids::Union{Vector{String},Nothing})

Compute global means for datasets in `data` with `ids`.
    
If `ids` is not specified, compute global means for all datasets in `data`.
"""
function addGlobalMeans!(
    data::DataMap; ids::Union{Vector{String}, Nothing}=nothing
 )
    if isnothing(ids)
        # NOTE: collect is important here, otherwise 'ids' changes when new 
        # keys are added in for-loop below to the data dictionary!!
        ids = collect(keys(data))
    elseif any(map(x -> !haskey(data, x), ids))
        throw(ArgumentError("addGlobalMeans!: There is data missing for some $ids !"))
    end
    for id in ids
        addDiagnostic!(data, computeGlobalMeans, id)
    end
    return nothing
end


function computeAnomalies(
    orig_data::YAXArray, ref_data::YAXArray; stats::String="ANOM"
)
    dimension, models = getDimsModel(orig_data)
    if !hasdim(ref_data, dimension) || models != Array(dims(ref_data, dimension))
        throw(ArgumentError("Original and reference data must contain exactly the same models!"))
    end
    if orig_data.properties["units"] != orig_data.properties["units"]
        @warn "Data and reference data are given in different units! NO ANOMALIES computed!"
        return nothing
    end
    # TODO: add check that datasets have the same units!
    anomalies_metadata = deepcopy(orig_data.properties)
    anomalies_metadata["_statistic"] = stats
    anomalies_id = buildMetaDataID(anomalies_metadata)
    anomalies_metadata["_id"] = anomalies_id
    anomalies_metadata["_ref_data_id"] = ref_data.properties["_id"]
    anomalies_metadata["_orig_data_id"] = orig_data.properties["_id"]
    
    anomalies = @d orig_data .- ref_data
    anomalies = YAXArray(
        dims(orig_data), Array(anomalies), anomalies_metadata
    )
    return anomalies
end


"""
    addAnomalies!(
        datamap::DataMap id_data::String, id_ref::String; stats::String="ANOM"
    )

Add entry to `datamap` with difference between `datamap` at `id_data` and `id_ref`.

The id of the original data and of the reference data is added to the metadata.
"""
function addAnomalies!(
    data::DataMap, id_data::String, id_ref::String; stats::String="ANOM"
)
    if !haskey(data, id_data)
        throw(ArgumentError("The given DataMap does not have key $id_data"))
    elseif !haskey(data, id_ref)
        throw(ArgumentError("The given DataMap does not have key $id_ref"))
    end
    addDiagnostic!(data, computeAnomalies, id_data, data[id_ref]; stats)
    return nothing
end


"""
    addAnomaliesGM!(data::DataMap, ids_data::Vector{String})

Compute anomalies of datasets in `data` with ids in `ids_data` with respect to 
their global means.
"""
function addAnomaliesGM!(data::DataMap, ids_data::Vector{String})
    # only compute global means for data for which it isnt already there
    gms_ids = similar(ids_data)
    for (i, id) in enumerate(ids_data)
        new_stats = hasdim(data[id], :time) ? "GM-ts" : "GM"
        gms_ids[i] = replace(id, data[id].properties["_statistic"] => new_stats)
    end
    ids = Vector{String}()
    for idx in range(1, length(gms_ids))
        if !haskey(data, gms_ids[idx])
            push!(ids, ids_data[idx])
        end
    end
    if !isempty(ids)
        addGlobalMeans!(data; ids)
    end
    @info "add anomalies wrt global mean..."
    for (i, id) in enumerate(ids_data)
        id_ref = gms_ids[i]
        addAnomalies!(
            data, id, id_ref; stats="ANOM-" * string(split(id_ref, "_")[2])
        )
    end
   return nothing
end


"""
    computeTempSTD(data::YAXArray, trend::YAXArray)

Compute the standard deviation of the temporal detrended timeseries `data` as 
STD^{t_2}_t=t_1 (X_l^t) - X_l^{TREND} where l is a (lon, lat)-position and 
TREND is the prediction of the ordinary least squares fit of the data (i.e. 
`trend` must have the same size as `data`).
"""
function computeTempSTD(data::YAXArray, trend::YAXArray)
    if otherdims(data, :time) != otherdims(trend, :time)
        msg = "Temporal Standard deviation can only be computed for timeseries data and trend data with identical dimensions!"
        throw(ArgumentError(msg))
    end
    meta_new = deepcopy(data.properties)
    meta_new["_id"] = replace(meta_new["_id"], meta_new["_statistic"] => "STD-temp")
    meta_new["_statistic"] = "STD-temp"

    if hasdim(trend, :time)
        diffs = @d data .- trend
        stds = dropdims(std(diffs, dims=:time), dims=:time)
    else
        # TODO
        @warn "computation of temporal trend with just slopes not yet implemented"
        throw(UndefVarError(:stds))
    end
    return YAXArray(otherdims(data, :time), stds.data, meta_new)
end


"""
    addTempSTD!(data::DataMap; statistic::String="CLIM-ann")

Compute temporal standard deviation (computeTempSTD) for every dataset with 
`statistic` in its id (default: CLIM-ann) and add the result to `data`. 

The id of computed temporal standard deviation is the same as before, but with 
`statistic` replaced by 'TREND'.
"""
function addTempSTD!(data::DataMap; statistic::String="CLIM-ann")
    ids = filter(id -> data[id].properties["_statistic"] == statistic, keys(data))
    if isempty(ids)
        @warn "No data for statistic $statistic found!"
        return nothing
    end
    for id in ids
        trend_id = replace(id, statistic=>"TREND-pred")
        get!(data, trend_id, computeLinearTrend(data[id]; full_predictions=true))
        addDiagnostic!(data, computeTempSTD, id, data[trend_id])
    end
    return nothing
end


# TODO: add fn to compute climatologies
# """
#     addClimatologies!(
#     datamap::DataMap, ids_data::Vector{String}, start_year, end_year
# )

# Add entry to 'datamap' with difference between 'datamap' at id_data and id_ref.

# The id of the original data and of the reference data is added to the metadata.
# """
# function addClimatologies!(data::DataMap, ids_data::Vector{String}, start_year, end_year)
#     anomalies = computeAnomalies(data[id_data], data[id_ref]; stats)
#     data[anomalies.properties["_id"]] = anomalies
#     return nothing
# end