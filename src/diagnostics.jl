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

    @timeit TO "build mask" begin 
        spatial = [1,2]
        other = filter(x -> !(x in spatial), 1:ndims(data))
        mask = falses(size(data)...)
        for I in CartesianIndices(axes(data)[other])
            indices = (:,:,I.I...)
            @views mask[indices...] .= ismissing.(@views data[indices...])
        end
    end
    aw_mat = areaWeightMatrix(latitudes, Array(mask)) # is normalized, lon x lat
    # make sure that, if provided, units are identical across models 
    units = get(data.properties, "units", nothing)
    if !isnothing(units) && length(unique(units)) != 1
        if any(x -> x != "degC" && x != "K", units)
            throw(ArgumentError("Data for all models must be defined in the same units to compute global mean! Only degC and K are converted into degC."))
        end
        data = kelvinToCelsius(data)
    end

    temp = aw_mat .* Array(data)
    gms = @timeit TO "last step" dropdims(mapslices(x -> sum(skipmissing(x)), temp; dims=(1, 2)); dims=(1,2))
    return YAXArray(otherdims(data, (:lon, :lat)), Array(gms), deepcopy(data.properties))
end

"""
# Arguments:
- `data::AbstractArray`: dimensions: lon, lat, model
"""
function globalMeans(data::AbstractArray, latitudes::AbstractArray)
    mask = Bool.(mapslices(model -> ismissing.(model), data; dims=(1,2)))
    aw_mat = areaWeightMatrix(latitudes, Array(mask)) # is normalized, lon x lat
    gms = mapslices(x -> sum(skipmissing(x)), aw_mat .* Array(data); dims=(1, 2))
    return dropdims(gms; dims=(1,2))
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
    units_orig = get(data.properties, "units", nothing)
    units_ref = get(ref_data.properties, "units", nothing) 
    if !(isnothing(units_orig) || isnothing(units_ref)) && units_ref != units_orig
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
    throwErrorIfDimMissing(data, dimension)
    standard_devs = dropdims(Statistics.std(data, dims = dimension), dims = dimension)
    return YAXArray(otherdims(data, dimension), standard_devs, deepcopy(data.properties))
end


function climatology(data::YAXArray)
    throwErrorIfDimMissing(data, :time)
    return dropdims(mean(data; dims=:time); dims=:time)
end



"""
    function gregoryECS(data_tas::YAXArray, data_rtmt::YAXArray)

Approximate a model's Equilibriate Climate Sensitivity (ECS) value using the Gregory method,
based on data from the 4xCO2 and piControl experiment.

# Arguments:
- `data_tas::YAXArray`: timeseries of global means of tas anomalies wrt piControl experiment
for 4xCO2 experiment. Must have dimension :model
- `data_rtmt::YAXArray`: timeseries of global means of rtmt data for 4xCO2 experiment. Must
have dimension :model
"""
function gregoryECS(data_tas::YAXArray, data_rtmt::YAXArray)
    models = dims(data_tas, :model)
    estimated_ecs = YAXArray(
        (models,),
        zeros(length(models))
    )
    for m in models
        lr_model = lr.linregress(data_tas[model=At(m)].data, data_rtmt[model=At(m)].data)
        slope, bias = lr.coef(lr_model)
        estimated_ecs[model = At(m)] = (-1) * bias / (2*slope)
    end
    return estimated_ecs
end