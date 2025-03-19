using DimensionalData
using GLM
using YAXArrays



"""
    getLinearTrend(data::YAXArray; full_predictions::Bool=true)

Compute linear trend as ordinary least squares for timeseries data.

If `data` has dimensions 'lon' and 'lat', linear trend is computed for each 
position separately.

# Arguments:
- `data::YAXArray`: must have dimension :time and one of :model, :member.
- `full_predictions::Bool`: if false, return only slope, otherwise return 
predicted value.
"""
function getLinearTrend(data::YAXArray; full_predictions::Bool=true)
    if !hasdim(data, :time) || (!hasdim(data, :model) && !hasdim(data, :member))
        msg = "Linear trend can only be computed for timeseries data (with dimension :time)!"
        throw(ArgumentError(msg))
    end
    has_grid = hasdim(data, :lon) && hasdim(data, :lat)
    x = Dates.year.(Array(data.time))
    meta = deepcopy(data.properties)
    stat = full_predictions ? "TREND-pred" : "TREND"
    meta["_id"] = replace(meta["_id"], meta["_statistic"] => stat)
    meta["_statistic"] = stat
    
    dimensions = full_predictions ? dims(data) : otherdims(data, :time)
    trends = YAXArray(dimensions, Array{eltype(data)}(undef, size(dimensions)), meta)

    is_level_model = hasdim(data, :model)
    models = is_level_model ? dims(data, :model) : dims(data, :member)
    for m in models
        if has_grid
            data_model = is_level_model ? data[model = At(m)] : data[member = At(m)]
            for (idx, y) in pairs(eachslice(data_model, dims=(:lon, :lat)))
                ols = lm(@formula(Y~X), DataFrame(X=x, Y=Float64.(Array(y))))
                lon_idx, lat_idx = Tuple(idx)
                if full_predictions
                    if is_level_model 
                        trends[lon=lon_idx, lat=lat_idx, model=At(m), time=:] .= predict(ols)
                    else 
                        trends[lon=lon_idx, lat=lat_idx, member=At(m), time=:] .= predict(ols)
                    end
                else
                    if is_level_model
                        trends[lon=lon_idx, lat=lat_idx, model=At(m)] = coef(ols)[2]
                    else
                        trends[lon=lon_idx, lat=lat_idx, member=At(m)] = coef(ols)[2]
                    end
                end
            end
        else
            y = is_level_model ? Array(data[model = At(m)]) : Array(data[member = At(m)])
            ols = lm(@formula(Y~X), DataFrame(X=x, Y=y))
            if full_predictions
                if is_level_model
                    trends[model = At(m)] .= predict(ols)
                else
                    trends[member = At(m)] .= predict(ols)
                end
            else
                if is_level_model
                    trends[model = At(m)] = coef(ols)[2]
                else 
                    trends[member = At(m)] = coef(ols)[2]
                end
            end
        end
    end
    return trends
end


"""
    addLinearTrend!(
        data::DataMap;     
        statistic::String="CLIM-ann", 
        ids_ts::Vector{String}=Vector{String}(), 
        full_predictions::Bool=true
    )

Add computed linear trend (getLinearTrend) for all datasets in `data` with an 
id in `ids_ts`, or if `ids_ts` is not given for all annual climatologies 
(stats=CLIM-ann) or with given `statistic` in `data`.
"""
function addLinearTrend!(
    data::DataMap; 
    statistic::String="CLIM-ann", 
    ids_ts::Vector{String}=Vector{String}(), 
    full_predictions::Bool=true
)
    if isempty(ids_ts)
        ids_ts = filter(id -> data[id].properties["_statistic"] == statistic, keys(data)) 
    end
    for id in ids_ts
        @info "add trend for $id"
        # TODO: there may be problems when data contains missing values!
        trend = getLinearTrend(data[id]; full_predictions)
        data[trend.properties["_id"]] = trend
    end
    return nothing
end



"""
    computeGlobalMeans(data::AbtractArray)

Return a YAXArray with area-weighted global means for each model in `data`. 

Missing data is accounted for in the area-weights. 

# Arguments:
- `data::YAXArray`: has at least dimensions 'lon', 'lat' and possibly 
'member' or 'model'.
"""
function computeGlobalMeans(data::YAXArray)
    longitudes = Array(dims(data, :lon))
    latitudes = Array(dims(data, :lat))
    masks = ismissing.(data)
    dimension = hasdim(data, :member) ? :member : hasdim(data, :model) ? :model : nothing
    meta = makeMetadataGMS(data.properties)
    if !isnothing(dimension)
        models = Array(dims(data, dimension))
        global_means = YAXArray(
            (Dim{dimension}(models),),
            Array{Union{Float64, Missing}}(undef, length(models)), 
            meta
        )
        for model in models
            mask = getAtModel(masks, dimension, model)
            global_mean = missing
            if any(mask .== false)
                area_weights = makeAreaWeightMatrix(longitudes, latitudes; mask)
                vals = dimension == :model ? data[model = At(model)] : data[member = At(model)]
                global_mean = Statistics.sum(skipmissing(vals .* area_weights))
            end
                putAtModel!(global_means, dimension, model, global_mean)
        end
    else 
        area_weights = makeAreaWeightMatrix(longitudes, latitudes; mask=masks)
        # TODO: check that where data is missing area_weights should be 0!
        global_means = Statistics.sum(skipmissing(data .* area_weights))
    end
    return global_means
end


"""
    addGlobalMeans!(data::DataMap; ids::Union{Vector{String},Nothing})

Compute global means for datasets in `data` with `ids`.
    
If `ids` is not specified, compute global means for all datasets in `data`.


# Arguments:
- `ids::Union{Vector{String},Nothing}`: ids for which global means are computed.
"""
function addGlobalMeans!(
    data::DataMap; ids::Union{Vector{String},Nothing}=nothing
 )
    if isnothing(ids)
        # NOTE: collect is important here, otherwise 'ids' changes when new 
        # keys are added in for-loop below to the data dictionary!!
        ids = collect(keys(data))
    elseif any(map(x -> !haskey(data, x), ids))
        throw(ArgumentError("addGlobalMeans!: There is data missing for some $ids !"))
    end
    for id in ids
        @info "add global means for $id"
        dat = copy(data[id])
        gms = hasdim(dat, :time) ? computeGlobalMeansTS(dat) : computeGlobalMeans(dat)
        data[gms.properties["_id"]] = gms
    end
    return nothing
end


function computeGlobalMeansTS(data::YAXArray)
    dimension = hasdim(data, :member) ? :member : hasdim(data, :model) ? :model : nothing
    nb_timesteps = length(dims(data, :time))
    meta = makeMetadataGMS(data.properties)
    global_means = !isnothing(dimension) ? 
        YAXArray(
            (dims(data, dimension), dims(data, :time)),
            Array{Union{Float64, Missing}}(undef, (length(dims(data, dimension)), nb_timesteps)),
            meta
        ) :
        YAXArray(
            dims(data, :time),
            Array{Union{Float64, Missing}}(undef, (nb_timesteps,)),
            meta
        )
    for t in dims(data, :time)
        global_means[time = At(t)] = computeGlobalMeans(data[time = At(t)])
    end
    return global_means
end


function computeAnomalies(orig_data::YAXArray, ref_data::YAXArray; stats::String="ANOM")
    dimension = hasdim(orig_data, :member) ? :member : :model
    models = dims(orig_data, dimension)
    if !hasdim(ref_data, dimension) || models != dims(ref_data, dimension)
        throw(ArgumentError("Original and reference data must contain exactly the same models!"))
    end
    is_timeseries = hasdim(orig_data, :time)

    anomalies_metadata = deepcopy(orig_data.properties)
    anomalies_metadata["_statistic"] = stats
    anomalies_id = buildMetaDataID(anomalies_metadata)
    anomalies_metadata["_id"] = anomalies_id
    anomalies_metadata["_ref_data_id"] = ref_data.properties["_id"]
    anomalies_metadata["_orig_data_id"] = orig_data.properties["_id"]

    arr = Array{Union{Missing, Float64}}(undef, size(dims(orig_data))...)
    anomalies = YAXArray(dims(orig_data), arr, anomalies_metadata)
    for m in models
        if is_timeseries
            times = dims(orig_data, :time)
            for t in times
                if dimension == :model
                    anomalies[time=At(t), model=At(m)] .= orig_data[time=At(t), model=At(m)] .- ref_data[time=At(t), model=At(m)]
                else
                    anomalies[time=At(t), member=At(m)] .= orig_data[time=At(t), member=At(m)] .- ref_data[time=At(t), member=At(m)]
                end
            end
        else
            if dimension == :model
                anomalies[model=At(m)] .= orig_data[model=At(m)] .- ref_data[model=At(m)]
            else
                anomalies[member=At(m)] .= orig_data[member=At(m)] .- ref_data[member=At(m)]
            end
        end
    end
    return anomalies
end

"""
    addAnomalies!(datamap::DataMap id_data::String, id_ref::String)

Add entry to 'datamap' with difference between 'datamap' at id_data and id_ref.

The id of the original data and of the reference data is added to the metadata.
"""
function addAnomalies!(data::DataMap, id_data::String, id_ref::String; stats::String="ANOM")
    if !haskey(data, id_data)
        throw(ArgumentError("The given DataMap does not have key $id_data"))
    elseif !haskey(data, id_ref)
        throw(ArgumentError("The given DataMap does not have key $id_ref"))
    end
    anomalies = computeAnomalies(data[id_data], data[id_ref]; stats)
    data[anomalies.properties["_id"]] = anomalies
    return nothing
end


function addAnomaliesGM!(data::DataMap, ids_data::Vector{String})
    # only compute global means for data for which it isnt already there
    gms_ids = map(id -> replace(id, data[id].properties["_statistic"] => "GM"), ids_data)
    ids = Vector{String}()
    for idx in range(1, length(gms_ids))
        if !haskey(data, gms_ids[idx])
            push!(ids, ids_data[idx])
        end
    end
    if !isempty(ids)
        addGlobalMeans!(data; ids)
    end
    for id in ids_data
        @info "add anomalies wrt global mean for $id"
        addAnomalies!(
            data, id, replace(id, data[id].properties["_statistic"] => "GM"); 
            stats="ANOM-GM"
        )
    end
   return nothing
end



"""
    computeTempSTD(data::YAXArray, trend::YAXArray)

Compute the standard deviation of the temporal detrended timeseries `data`.
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
        stds = dropdims(std(data .- trend, dims=:time), dims=:time)
    else
        # TODO
        @warn "computation of temporal trend with just slopes not yet implemented"
        throw(UndefVarError(:stds))
    end
    return YAXArray(otherdims(data, :time), stds.data, meta_new)
end


function addTempSTD!(data::DataMap; statistic::String="CLIM-ann")
    ids = filter(id -> data[id].properties["_statistic"] == statistic, keys(data))
    for id in ids
        trend_id = replace(id, statistic=>"TREND-pred")
        @info "add temp std for data: $id and trend $trend_id"
        get!(data, trend_id, getLinearTrend(data[id]))
        standard_dev = computeTempSTD(data[id], data[trend_id])
        data[standard_dev.properties["_id"]] = standard_dev
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