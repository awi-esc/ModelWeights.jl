using DimensionalData
using GLM
using YAXArrays


"""
    computeLinearTrend(data::YAXArray; full_predictions::Bool=true)

Compute linear trend as ordinary least squares for timeseries data.

In the metadata of the returned array '-TREND' (if `full_predictions`=false) or 
'-TREND-pred' (if `full_predictions`=true) is added to the original metadata 
entry of `data` for '_statistic' and the new metadata entry for '_id' is adapted accordingly.

# Arguments:
- `data::YAXArray`: must have dimension :time.
- `full_predictions::Bool`: if false, return only slope, otherwise return predicted value.
"""
function computeLinearTrend(
    data::YAXArray; stats_name::String="", full_predictions::Bool = false, only_newid::Bool=false
)
    if !hasdim(data, :time)
        msg = "Linear trend can only be computed for timeseries data (with dimension :time)!"
        throw(ArgumentError(msg))
    end
    if isempty(stats_name)
        stats_name = full_predictions ? "TREND-pred" : "TREND"
    end
    x = Dates.year.(Array(data.time))
    meta = deepcopy(data.properties)
    meta["_statistic"] = stats_name
    meta["_id"] = buildMetaDataID(meta)
    if only_newid
        return meta["_id"]
    else
        function fn(y)
            y = Array(y)
            indices = findall(.!(ismissing.(y)))
            # TODO: this is slow and should be changed!
            ols = lm(@formula(Y ~ X), DataFrame(X = x[indices], Y = Float64.(y[indices])))
            if full_predictions
                y_hat = Array{Union{Missing,Float64}}(undef, size(y))
                y_hat[indices] = predict(ols)
            else
                y_hat = coef(ols)[2]
            end
            return y_hat
        end
        trends = mapslices(fn, data; dims = (:time,))
        return YAXArray(dims(trends), Array(trends), meta)
    end
end


function detrendTimeseries(data::YAXArray; stats_name::String="detrended-ts", only_newid=false)
    meta_new = deepcopy(data.properties)
    meta_new["_statistic"] = stats_name
    meta_new["_id"] = buildMetaDataID(meta_new)
    if only_newid
        return meta_new["_id"]
    else
        trend = computeLinearTrend(data; full_predictions = true)
        diffs = @d data .- trend
    end
    return YAXArray(dims(data), diffs, meta_new)
end