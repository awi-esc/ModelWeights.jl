module Timeseries

using Dates
using DimensionalData
using LinearRegression
using Missings
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
    trend = linearTrend(data; full_predictions = true)
    return detrend(data, trend)
end

function detrend(data::YAXArray, trend::YAXArray)
    data_mat = YAXArray(data.axes, Array(data), deepcopy(data.properties))
    trend_mat = YAXArray(trend.axes, Array(trend), deepcopy(trend.properties))
    return @d data_mat .- trend_mat
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
            Data.subsetMeta!(df.properties, indices_keep)
        end
    end
    df.properties["timerange"] = join(string.([start_year, end_year]), "-")
    warnIfOutOfTimerange(df, start_year, end_year)
    Data.warnIfhasMissing(df; name = "timeseries data")
    return df
end

"""
    function filterTimeseries(data::YAXArray, indices::AbstractArray{<:Int})

# Arguments:
- `data::YAXArray`:
- `indices::AbstractArray{<:Int}`: start and end indices
"""
function filterTimeseries(
    data::YAXArray,
    indices::AbstractArray{<:Int}; 
    relative::Bool=false
)
    Data.throwErrorIfDimMissing(data, :time)
    times = Array(data.time)
    if relative
        indices_ts = get(data.properties, "timeseries-start-end", nothing)
        if isnothing(indices_ts) 
            indices_ts = addStartEnd!(data)
        end
        ts_lengths = indices_ts[:,2] .- indices_ts[:,1] .+ 1
        ts_lengths[ts_lengths .== 1] .= 0

        n = indices[2]-indices[1] + 1
        dim_time = map(x -> DateTime(x), 1 : n)
        dimensions = Vector(undef, ndims(data))
        i_time = Data.indexDim(data, :time)
        i_model = Data.indexDim(data, :model)
        for (i,d) in enumerate(dims(data))
            dimensions[i] = i == i_time ? Dim{:time}(dim_time) : d
        end
        s = collect(size(data))
        s[i_time] = n
        s_temp = copy(s)
        dat = YAXArray(Tuple(dimensions), allowmissing(zeros(Tuple(s))))
        for (i, m) in enumerate(data.model)
            if ts_lengths[i] == 0
                s_temp[i_model] = 1
                vals = fill(missing, s_temp...)
            else
                ts_start = indices_ts[i,1] + indices[1] - 1
                ts_end = ts_start + n - 1
                #TODO: this fils if ts_start or ts_end exceed boundaries
                vals = data[model = At(m), time = ts_start : ts_end]
            end
            dat[model = At(m)] = vals
        end
    else
        dat = data[time = indices[1] : indices[2]]
        # TODO: add case when indices exceeds range
        dat = YAXArray(dims(dat), Array(dat), deepcopy(dat.properties))
        start_year = times[1]
        end_year = times[end]
        dat.properties["timerange"] = join(string.([start_year, end_year]), "-")
    end
    return dat
end


function addStartEnd!(data::YAXArray)
    i_time = Data.indexDim(data, :time)
    i_model = Data.indexDim(data, Data.modelDim(data))
    models = lookup(data, i_model)
    indices = zeros(Int, length(models), 2)
    times = data.time
    n_timesteps = length(times)
    nb_dims_larger_2 = ndims(data) > 2

    model_dim_last = i_model > i_time

    arr = data.data
    for (i, model) in enumerate(models)
        ts = 1;
        i_start = 0; i_end = 0;
        while ts < n_timesteps && i_start == 0
            #slice = @view data[model = At(model), time=ts]
            slice = model_dim_last ? selectdim(arr, i_model, i) : selectdim(arr, i_time, ts)
            slice = model_dim_last ? selectdim(slice, i_time, ts) : selectdim(slice, i_model, i)
            all_missing = nb_dims_larger_2 ? all(ismissing, slice) : ismissing(slice.data[])
            if all_missing
                ts += 1
            else
                i_start = ts
                ts += 1
            end
        end

        if ts == n_timesteps
            # all missing
        else
            while ts < n_timesteps && i_end == 0
                #slice = @view data[model = At(model), time=ts]
                slice = model_dim_last ? selectdim(arr, i_model, i) : selectdim(arr, i_time, ts)
                slice = model_dim_last ? selectdim(slice, i_time, ts) : selectdim(slice, i_model, i)

                all_missing = nb_dims_larger_2 ? all(ismissing, slice) : ismissing(slice.data[])
                if all_missing
                    i_end = ts-1
                else
                    ts += 1
                end
            end
            if ts == n_timesteps
                i_end = n_timesteps
            end
        end
        indices[i, 1] = i_start
        indices[i, 2] = i_end
    end
    data.properties["timeseries-start-end"] = indices
    return indices
end



end # end module