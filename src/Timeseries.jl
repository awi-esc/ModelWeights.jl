module Timeseries

using DataFrames
using Dates
using DimensionalData
using GLM
using YAXArrays

using ..Data

"""
    linearTrend(data::YAXArray; full_predictions::Bool=true)

Compute linear trend as ordinary least squares for timeseries data.

In the metadata of the returned array '-TREND' (if `full_predictions`=false) or 
'-TREND-pred' (if `full_predictions`=true) is added to the original metadata 
entry of `data` for 'statistic' and the new metadata entry for 'id' is adapted accordingly.

# Arguments:
- `data::YAXArray`: must have dimension :time.
- `full_predictions::Bool`: if false, return only slope, otherwise return predicted value.
"""
function linearTrend(data::YAXArray; full_predictions::Bool = false)
    if !hasdim(data, :time)
        msg = "Linear trend can only be computed for timeseries data (with dimension :time)!"
        throw(ArgumentError(msg))
    end
    x = Dates.year.(Array(data.time))
    meta = deepcopy(data.properties)
    meta["statistic"] = full_predictions ? "TREND-pred" : "TREND"
    meta["id"] = Data.buildMetaDataID(meta)
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


function detrend(data::YAXArray)
    meta_new = deepcopy(data.properties)
    meta_new["statistic"] = join([meta_new["statistic"], "detrended-ts"], "-")
    meta_new["id"] = Data.buildMetaDataID(meta_new)
    trend = linearTrend(data; full_predictions = true)
    diffs = @d data .- trend
    return YAXArray(dims(data), diffs, meta_new)
end


function warnIfOutOfTimerange(df::YAXArray, start_year::Int, end_year::Int)
    timesteps = dims(df, :time)
    data_start = Dates.year(timesteps[1])
    if data_start != start_year
        @warn "start_year for is : $(data_start) instead of $start_year !"
    end
    data_end = Dates.year(timesteps[end])
    if data_end != end_year
        @warn "end_year is : $(data_end) instead of $end_year"
    end
    return nothing
end


"""
    function filterTimeseries(
        data::YAXArray,
        start_year::Number,
        end_year::Number;
        only_models_non_missing_vals::Bool = true
    )

# Arguments:
- `data::YAXArray`:
- `start_year::Number`: set to -Inf to not filter for start time.
- `end_year::Number`: set to Inf to not filter for end time.
- `only_models_non_missing_vals::Bool = true`:
"""
function filterTimeseries(
    data::YAXArray,
    start_year::Number,
    end_year::Number;
    only_models_non_missing_vals::Bool = true
)
    Data.throwErrorIfDimMissing(data, :time)
    try
        df = data[time = Where(
            x -> Dates.year(x) >= start_year && Dates.year(x) <= end_year,
        )]
        df = YAXArray(dims(df), Array(df), deepcopy(df.properties))
    catch
        @warn "No data found in between $start_year and $end_year !"
        return nothing
    end
    if only_models_non_missing_vals
        dim_symbol = Data.modelDim(df)
        models_missing_vals = mapslices(x -> any(ismissing.(x)), df, dims = otherdims(df, dim_symbol))
        n_models = length(dims(models_missing_vals, dim_symbol))
        indices_missing = findall(x -> x == true, models_missing_vals)
        if !isempty(indices_missing)
            models_missing = Data.getByIdxModel(models_missing_vals, dim_symbol, indices_missing)
            models_missing = dims(models_missing, dim_symbol).val
            df = dim_symbol == :model ? df[model = Where(x -> !(x in models_missing))] :
                df[member = Where(x -> !(x in models_missing))]
            # update metadata too
            indices_keep = filter(x -> !(x in indices_missing), 1:n_models)
            Data.summarizeMeta!(df.properties, indices_keep)
        end
    end
    df.properties["timerange"] = join(string.([start_year, end_year]), "-")
    warnIfOutOfTimerange(df, start_year, end_year)
    Data.warnIfhasMissing(df; name = "timeseries data")
    return df
end


"""
    addFilteredTimeseries!(
        data::DataMap, 
        start_year::Number, 
        end_year::Number;
        ids_ts::Vector{String}=Vector{String}(),
        new_alias::String,
        only_models_non_missing_vals::Bool = true
    )
"""
function addFilteredTimeseries!(
    data::DataMap,
    start_year::Number,
    end_year::Number;
    ids_ts::Vector{String} = Vector{String}(),
    new_alias::String = "",
    only_models_non_missing_vals::Bool = true,
)
    if isempty(ids_ts)
        @info "filter all datasets with :time dimension from $start_year to $end_year"
        ids_ts = filter(id -> hasdim(data[id], :time), collect(keys(data)))
    else
        indices_keys = map(id -> haskey(data, id), ids_ts)
        ids_absent = ids_ts[.!(indices_keys)]
        if !isempty(ids_absent)
            if length(ids_absent) == length(ids_ts)
                throw(ArgumentError("No data found for ids: $(ids_ts)"))
            else
                @warn "No data found for ids" ids_absent
            end
        end
        ids_ts = ids_ts[indices_keys]
    end
    for id in ids_ts
        df = filterTimeseries(data[id], start_year, end_year; new_alias, only_models_non_missing_vals)
        data[df.properties["id"]] = df
    end
    return nothing
end


function addFilteredTimeseries!(
    data::ClimateData,
    start_year::Number,
    end_year::Number;
    ids_ts::Vector{String} = Vector{String}(),
    new_alias::String="",
    only_models_non_missing_vals::Bool = true,
)   
    @info "addFilteredTimeseries for Model data..."
    addFilteredTimeseries!(data.models, start_year, end_year; ids_ts, new_alias, only_models_non_missing_vals)
    @info "addFilteredTimeseries for Observational data..."
    addFilteredTimeseries!(data.obs, start_year, end_year; ids_ts, new_alias, only_models_non_missing_vals)
    return nothing
end

end # end module