module Timeseries

using Dates
using DimensionalData
using LinearRegression
using YAXArrays

using ..Data

"""
    linearTrend(data::YAXArray; full_predictions::Bool=true)

Compute linear trend as ordinary least squares for timeseries `data` and return slope 
(default) or predicted y for every element in `x` (if `full_predictions`=true)

# Arguments:
- `data::YAXArray`: must have dimension :time.
- `full_predictions::Bool`:
"""
function linearTrend(data::YAXArray; full_predictions::Bool = false)
    Data.throwErrorIfDimMissing(data, :time)
    timesteps = Dates.year.(Array(data.time))
    trends = mapslices(y -> computeLinearTrend(timesteps, Array(y); full_predictions), data; dims = (:time,))
    return full_predictions ? trends : dropdims(trends, dims=:OutAxis1)
end


"""
    computeLinearTrend(x::Vector, y::Vector; full_predictions::Bool = false)

Regress y on x and return slope (default) or predicted y for every element in `x` (if `full_predictions`=true)

# Arguments:
- `x::Vector`:
- `y:Vector`:
- `full_predictions::Bool`: if false, return slope, otherwise return all values predicted for `x`.
"""
function computeLinearTrend(x::Vector, y::Vector; full_predictions::Bool = false)
    indices = findall(.!(ismissing.(y)))
    lr = linregress(x[indices], y[indices])
    slope = LinearRegression.slope(lr)
    return full_predictions ? LinearRegression.bias(lr) .+ (slope .* x) : slope
end


function detrend(data::YAXArray)
    return @d data .- linearTrend(data; full_predictions = true)
end

function detrend(data::YAXArray, trend::YAXArray)
    return @d data .- trend
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