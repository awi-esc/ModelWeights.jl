"""
    nbModelMembers(data::YAXArray)

Return a dictionary mapping from models in `data` to the number of its members.
"""
function nbModelMembers(data::YAXArray)
    throwErrorIfDimMissing(data, :member)
    data = setLookupsFromMemberToModel(data, ["member"])
    models = collect(lookup(data, :model))
    return countMap(models)
end

"""
    convertToYAX(dm::DataMap; dim_name::Symbol = :diagnostic)

Convert a DataMap into a YAXArray. All entries in `dm` must have the same dimensions.
"""
function convertToYAX(dm::DataMap; dim_name::Symbol = :diagnostic)
    # TODO: add dimension checks!
    diagnostics = collect(keys(dm))
    dimensions = dims(dm[diagnostics[1]])

    data = YAXArray(
        (dimensions..., Dim{dim_name}(diagnostics)), 
        rand(size(dimensions)..., length(diagnostics))
    )
    for d in diagnostics        
        data[diagnostic = At(d)] .= dm[d]
    end
    return data
end


"""
    getAtModel(data::YAXArray, dimension::Symbol, model::String)

Return `data` where `dimension` (member or model) has value `model`.
"""
function getAtModel(data::YAXArray, dimension::Symbol, model::String)
    throwErrorIfDimMissing(data, dimension)
    return dimension == :model ? data[model=At(model)] : data[member=At(model)]
end


"""
    getByIdxModel(data::YAXArray, dimension::Symbol, indices::Vector)

Return `data` where `dimension` (member or model) has value `model`.
"""
function getByIdxModel(data::YAXArray, dimension::Symbol, indices::Vector)
    throwErrorIfDimMissing(data, dimension)
    return dimension == :model ? data[model=indices] : data[member=indices]
end


"""
    indexModel(data::YAXArray, model_dims::Tuple{Symbol}, indices::Vector{Int})

Return `data` at model dimensions `model_dims` at `indices`.
"""
function indexModel(
    data::YAXArray, model_dims::NTuple{N, Symbol}, indices::Vector{Int}
) where {N}
    data = YAXArray(data.axes, data.data, deepcopy(data.properties))
    dim_names = dimNames(data)
    index_vec = map(dim_names) do name 
        name in model_dims ? indices : Colon()
    end
    data = data[index_vec...] # this indexing materializes the data, i.e. its not lazy anymore!
    subsetMeta!(data.properties, indices)
    return data
end



""" 
    putAtModel!(data::YAXArray, dimension::Symbol, model::String, input)
"""
function putAtModel!(data::YAXArray, dimension::Symbol, model::String, input)
    throwErrorIfDimMissing(data, dimension)
    if dimension == :model
        data[model = At(model)] = input
    else
        data[member = At(model)] = input
    end
    return nothing
end


"""
    mergeMetaDataFromMultipleFiles(data::Vector{<:YAXArray})

Combine arrays in `data` into a single YAXArray with meta combined from all datasets into 
lists, with missing if key wasnt in a dataset.
"""
function mergeMetaDataFromMultipleFiles(data::AbstractVector{<:AbstractArray})
    n_files = length(data)
    meta_keys = unique(vcat(map(x -> collect(keys(x.properties)), data)...))
    meta_dict = Dict{String, Any}()
    for (i, ds) in enumerate(data)
        for key in meta_keys
            values = get!(meta_dict, key, repeat(Any[missing], outer=n_files))
            values[i] = get(ds.properties, key, missing)
        end
    end
    return meta_dict
end


"""
    combineModelsFromMultipleFiles(
        data::Vector{AbstractArray}; 
        model_names::Vector{String} = Vector{String}(),
        meta::Union{Dict{String, T}, Nothing} = nothing,
        dtype::String = "undef",
        new_dim::Symbol = :model,
        sorted::Bool = true
) where T <: Any

Combine `data` from different files into a single YAXArray. 
The meta data of the returned YAXArray is a combination of all values from all datasets 
using vectors with missing values for datasets that didn't have the respective metadata 
entry in their metadata, further combined with the properties in `meta` if provided.

All elements in `data` must share all dimensions except for time (if present).
For timeseries data, the time dimension may cover different ranges. In that case, the 
maximal possible timerange is used where missing values are added if the timepoint had not been 
defined for a model.

The combined YAXArray has the additional dimension `new_dim` (default: :model) with `names`
as values. If `names` is not provided, default values 'model1', 'model2, etc. are used 
instead.

If `sorted` is true, model dimension of returned data is sorted alphabetically and the 
vector entries in the metadata dictionary of the returned array are also sorted accordingly.
"""
function combineModelsFromMultipleFiles(
    data::AbstractVector{<:YAXArray}; 
    model_names::Vector{String} = Vector{String}(),
    meta::Union{Dict{String, T}, Nothing} = nothing,
    new_dim::Symbol = :model,
    sorted::Bool = true
) where T <: Any
    if isempty(data)
        @warn "Data vector is empty!"
        return nothing
    end
    data_sizes = unique(map(size, data))
    if length(data_sizes) != 1
        # dimensions must be identical except for time
        if !all(map(x -> hasdim(x, :time), data))
            msg = "Data does not have the same size across all models: $(data_sizes)"
            throw(ArgumentError(msg))
        else
            # if difference only in time, use maximal possible timeseries and add missing values
            alignTimeseries!(data)
            #overlapTimeseries!(data)
        end
    end
    use_default_names = isempty(model_names)
    if use_default_names
        model_names = map(x -> string(new_dim) * string(x), 1:length(data))
    end
    if sorted && !use_default_names
        sort_indices = sortperm(model_names)
        model_names = model_names[sort_indices]
        data = data[sort_indices]
    end
    meta_dict = mergeMetaDataFromMultipleFiles(data)    
    if !isnothing(meta)
        meta_dict["_meta"] = meta
    end
    # concatenatecubes doesnt work with 0-dimensional arrays...
    has_0_dims = isempty(dims(data[1]))
    if has_0_dims
        data = map(x -> YAXArray((Dim{:temp}([1]),), vec(x)), data)
    end
    dimData = concatenatecubes(data, Dim{new_dim}(model_names))
    if has_0_dims
        dimData = dimData[temp = 1]
    end
    dimData = YAXArray(dimData.axes, dimData, meta_dict)
    if length(model_names) != length(unique(model_names))
        duplicates = unique([m for m in model_names if sum(model_names .== m) > 1])
        # handle duplicates
        indices_remove = []
        models = collect(lookup(dimData, new_dim))
        for m in duplicates
            indices = findall(x -> x == m, models)
            push!(indices_remove, indices[2:end]...)
        end
        indices_keep = filter(x -> !(x in indices_remove), collect(1:length(models)))
        # TODO: fix to do this also with arbitrary new_dim name!
        dimData = getByIdxModel(dimData, new_dim, indices_keep)
    end
    return dimData
end


"""
    memberIDsFromPaths(all_paths::Vector{String})

For every path in `all_paths` return a string of the form modelname#memberID[_grid]
that identifies the corresponding model member. Filenames must follow the CMIP-standard
(proj_name_mip_exp_id_variable[_grid].nc).

For CMIP6 models the abbreviation of the grid is added to the model name.
"""
function memberIDsFromPaths(all_paths::Vector{String})
    all_filenames = split.(basename.(all_paths), "_")
    all_members = Vector{String}(undef, length(all_filenames))
    for (i, fn_parts) in enumerate(all_filenames)
        model = join(fn_parts[[2, 5]], MODEL_MEMBER_DELIM)
        # add grid to model name for CMIP6 models:
        if fn_parts[1] != "CMIP5"
            model = splitext(model * "_" * fn_parts[7])[1]
        end
        all_members[i] = model
    end
    return unique(all_members)
end


"""
    memberIDFromFilenameMeta(fn_meta::FilenameMeta)

Return a string of the form modelname#memberID[_grid] that identifies the corresponding model member.

For CMIP6 models the abbreviation of the grid is added to the model name.
"""
function memberIDFromFilenameMeta(fn_meta::FilenameMeta, mip_era::String)
    id = join([fn_meta.model, fn_meta.variant], MODEL_MEMBER_DELIM)
    if lowercase(mip_era) == "cmip6"
        id = join([id, fn_meta.grid], "_")
    end
    return id
end


"""
    getCMIPModelsKey(meta::Dict)

Return the respective key to retrieve model names in CMIP6 ('source_id') and 
CMIP5 ('model_id') data.

If both keys are present, 'source_id' used in CMIP6 models is returned, if none 
is present, throw ArgumentError.
"""
function getCMIPModelsKey(meta::Dict)
    attributes = keys(meta)
    if "source_id" in attributes
        if "model_id" in attributes
            msg1 = "Dictionary contains keys source_id and model_id, source_id is used! "
            @debug msg1 meta["source_id"] meta["model_id"]
        end
        return "source_id"
    elseif "model_id" in attributes
        return "model_id"
    else
        msg = "For CMIP data metadata must contain one of 'source_id' (pointing to names of CMIP6 models) or 'model_id' (CMIP5). This might happen if you try to load observational and model data from the same directory."
        throw(ArgumentError(msg))
    end
end


"""
    filterPathsSharedModels(
        paths::Vector{String}, 
        shared_models::Vector{String}, 
        fn_format::Union{Symbol, String}
    )

Every vector of paths in `all_paths` is filtered s.t. it only contains models or model 
members given in `shared_models`.

# Arguments:
- `paths`: contains paths to data files
"""
function filterPathsSharedModels(
    paths::Vector{String}, shared_models::Vector{String}, fn_format::Union{Symbol, String}
)
    if isempty(shared_models)
        @warn "No models shared across data!"
        return Vector{String}()
    end
    constraint = Dict("models" => shared_models)
    mask = maskFileConstraints(paths, fn_format, constraint)
    return paths[mask]
end

"""
    filterPathsSharedModels(
        all_paths::Vector{Vector{String}}, 
        level_shared::Level, 
        fn_format::Union{Symbol, String}
    )

# Arguments:
- `all_paths`: every entry refers to the paths to data files for the respective dataset
"""
function filterPathsSharedModels(
    all_paths::Vector{Vector{String}}, 
    level_shared::Level,
    fn_format::Union{Symbol, String}
)
    shared = sharedModels(all_paths, level_shared, fn_format)
    return map(paths -> filterPathsSharedModels(paths, shared, fn_format), all_paths)
end

function filterPathsSharedModels(
    all_paths::Vector{Vector{String}},
    level_shared::Union{String, Symbol}, 
    fn_format::Union{Symbol, String}
)
    return filterPathsSharedModels(all_paths, getLevel(level_shared), fn_format)
end


function modelsFromMemberIDs(members::AbstractVector{<:AbstractString}; uniq::Bool=false)
    models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
    return uniq ? unique(models) : models
end

function modelsFromMemberIDs(data::YAXArray; uniq::Bool=false)
    throwErrorIfDimMissing(data, :member)
    return modelsFromMemberIDs(string.(dims(data, :member)); uniq)
end


"""
    areaWeightMatrix(
        latitudes::AbstractVector{<:Number}, mask::YAXArray{T}
    ) where {T <: Union{Missing, Bool}}

Return matrix of size length(longitudes) x length(latitudes) with normalized area weights 
(approximated based on latitudes) and set to 0 at positions where mask is true.

# Arguments:
- `latitudes::AbstractVector{<:Number}`: to approximate area weights with cosine of latitudes.
- `mask::AbstractArray{T}`: lon, lat as first and second dimension respectively where lat must 
have same length as `latitudes`.
"""
function areaWeightMatrix(
    latitudes::AbstractVector{T}, mask::AbstractArray{Bool}
) where T<:Number #where {T <: Union{Missing, Bool}}
    if size(mask, 2) != length(latitudes)
        throw(ArgumentError("Second dimension of mask must be lat and equal $(length(latitudes)), found:$(size(mask))"))
    end
    area_weights = cos.(deg2rad.(latitudes))
    s = size(mask)
    n_lon = s[1]
    aw_mat = length(s) > 2 ? repeat(area_weights', n_lon, 1, s[3:end]...) : repeat(area_weights', n_lon)
    aw_mat[mask] .= 0
    return aw_mat ./ sum(aw_mat; dims=(1,2))
end


function checkDataStructure(path_data::String, dir_per_var::Bool)
    if !dir_per_var
        subdirs = filter(x -> isdir(joinpath(path_data, x)), readdir(path_data))
        if !("preproc" in subdirs)
            throw(ArgumentError("If dir_per_var is false, path_data must contain a directory named 'preproc'"))
        end
    end
    return nothing
end


"""
    checkConstraint(constraint::Union{Dict, Nothing})

Throw error if entries (except for level_shared) do not map to a vector of Strings.

Used to provide specific error messages of values that were not specified correctly.
"""
function checkConstraint(constraint::Union{Dict, Nothing})
    if !isnothing(constraint)
        for (k, v) in constraint
            if !(k in ["level_shared", "timeseries"]) && !isa(v, Vector{String})
                throw(ArgumentError("Constraint must map to a vector of Strings for key $k, found: $(typeof(v))"))
            end
            if k == "timeseries" && !isa(v, Dict)
                throw(ArgumentError("Constraint timeseries must map to a Dictionary!"))
            end
        end
    end
    return nothing
end


"""
    joinDicts(v::Dict...; warn_msg::String="")
"""
function joinDicts(v::Dict...; warn_msg::String="")
    result = typeof(v[1])()
    for dm in v
        shared_keys = sharedKeys(result, dm)
        if length(shared_keys) > 0
            dups = filter(!isempty, map(shared_keys) do x
                    dm[x] != result[x] ? (x, result[x], dm[x]) : ()
                end
            )
            if !isempty(dups)
                @warn warn_msg dups
            end
        end
        result = merge(result, dm)
    end
    return result
end


function physicsFromMember(member::String)
    regex = r"(p\d+)(f\d+)?(_.*)?$"
    return String.(match(regex, member).captures[1])
end


function joinDataMaps(v::Union{Dict{String, MetaData}, DataMap}...; warn_msg::String="")
    return joinDicts(v...; warn_msg)
end


"""
    writeDataToDisk(data, target_path::String)

Save `data` as Julia obj if `target_path` has ending '.jld2', otherwise save as binary.
If file at `target_path` already exists, timestamp is added if `overwrite` is false (default).
"""
function writeDataToDisk(data, target_path::String; overwrite::Bool = false)
    target_path = overwrite ? target_path : Data.individuatePath(target_path)
    if endswith(target_path, ".jld2")
        jldsave(target_path; data = data)
    else
        serialize(target_path, data)
    end
    @info "saved data to: $(target_path)"
    return nothing
end


"""
    readDataFromDisk(target_path::String; variable::String="")

Load data from `target_path`. If `target_path` ends with '.jld2', `variable` 
must be specified, otherwise  data is assumed to be binary.
"""
function readDataFromDisk(target_path::String; variable::String = "data")
    if !isfile(target_path)
        throw(ArgumentError("There does not exist a file at: $(target_path)"))
    end
    if endswith(target_path, ".jld2")
        if isempty(variable)
            throw(ArgumentError("To load .jld2 data, specify argument variable!"))
        end
        f = jldopen(target_path, "r")
        data = f[variable]
        close(f)
    else
        data = deserialize(target_path)
    end
    return data
end


"""
    setLookupsFromMemberToModel(data::YAXArray, dim_names::Vector{String})

Change the lookup values for the dimension 'member' to refer to the models, i.e.
they are not unique anymore. This is done in preparation to group the data by
the different models.

# Arguments:
- `data::YAXArray`: has at least dimensions in `dim_names`.
- `dim_names::Vector{String}`: names of dimensions to be changed, e.g. 'member', 
'member1' (would be changed to 'model', 'model1').
"""
function setLookupsFromMemberToModel(data::YAXArray, dim_names::Vector{String})
    # TODO: add checks for dim names!
    n_dims = length(dim_names)
    for (i, dim) in enumerate(dim_names)
        unique_members = Array(dims(data, Symbol(dim)))
        models = map(x -> String.(split(x, MODEL_MEMBER_DELIM)[1]), unique_members)
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model"
        data = setDim(data, dim, new_dim_name, models)
    end
    return data
end


"""
    membersToModels(members::AbstractArray{String})

# Arguments:
- `members::AbstractArray{String}`: unique names of model members with model name followed by separator followed by unique identifiers.
"""
function membersToModels(members::AbstractArray{String})
    return map(x -> String.(split(x, MODEL_MEMBER_DELIM)[1]), members)
end


"""
    subsetMeta(meta::Dict, indices::Vector{<:Int}; simplify::Bool = false)

If simplify is true, use single value when all remaining elements in a vector in `meta` are 
identical, otherwise just keep vector.

# Arguments:
- `meta`:
- `indices`: indices to remain in vectors mapped to in `meta`. 
"""
function subsetMeta!(meta::Dict, indices::Vector{<:Int}; simplify::Bool = false)    
    for (k, v) in meta
        if isa(v, Vector)
            vals = v[indices]
            if simplify
                vals = unique(v[indices])
                vals = length(vals) == 1 ? vals[1] : vals
            end
        else 
            vals = deepcopy(v)
        end
        meta[k] = vals
    end
    return nothing
end


"""
    subsetMeta(meta::Dict, indices::Vector{<:Int}; simplify::Bool = false)

If simplify is true, use single value when all remaining elements in a vector in `meta` are 
identical, otherwise just keep vector.

# Arguments:
- `meta`:
- `indices`: indices to remain in vectors mapped to in `meta`. 
"""
function subsetMeta(meta::Dict, indices::Vector{<:Int}; simplify::Bool = false)    
    meta_new = Dict()    
    for (k, v) in meta
        if isa(v, Vector)
            vals = v[indices]
            if simplify
                vals = unique(v[indices])
                vals = length(vals) == 1 ? vals[1] : vals
            end
        else 
            vals = deepcopy(v)
        end
        meta_new[k] = vals
    end
    return meta_new
end

"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Summarize vectors in `meta`, refering to different models (members), such that 
each vector only contains N entries (N=number of models (i.e. without unique members)).

If the metadata for members of a model differ across members, the respective
entry in the vector will be a vector itself.
"""
function updateGroupedDataMetadata(
    meta::Dict, grouped_data::DimensionalData.DimGroupByArray,
)
    meta_new = filter(((k, v),) -> k == "member_names" || !(v isa Vector), meta)
    attributes = filter(x -> meta[x] isa Vector && x != "member_names", keys(meta))
    attribs_diff_across_members = []
    # iterate over attributes that are vectors, thus different for the different
    # members or models
    for key in attributes
        for (i, model) in enumerate(dims(grouped_data, :model))
            indices = findall(x -> x == model, meta["model_names"])
            vals = get!(meta_new, key, [])
            val_model = meta[key][indices]
            if length(unique(val_model)) != 1
                push!(vals, val_model)
                push!(attribs_diff_across_members, key)
            else
                # members of current ensemble all share the same value
                push!(vals, val_model[1])
            end
        end
    end
    if !isempty(attribs_diff_across_members)
        # TODO ignore those that are defenitely expected to differ
        @debug "metadata attributes that differ across model members (ok for some!)" unique(
            attribs_diff_across_members,
        )
    end
    return meta_new
end


"""
    metaDataFromESMValToolRecipes(
        base_path_configs::String; constraint::Union{Dict, Nothing} = nothing
    )

Read variable, statistic, experiment and timerange/alias values from ESMValTool recipes 
stored at `base_path_configs` into a vector of Dictionaries storing the respective readoff 
values.
"""
function metaDataFromESMValToolRecipes(
    base_path_configs::String;
    constraint::Union{Dict, Nothing} = nothing 
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join = true),
    )
    meta_attribs = Vector{MetaData}()
    for path_config in paths_configs
        config = YAML.load_file(path_config)
        data_all = config["diagnostics"]
        aliases = keys(data_all)

        for alias in aliases
            data = data_all[alias]["variables"]
            for (k, v) in data      
                variable = get(v, "short_name", "")
                var_long = get(v, "variable_long_name", "")
                preprocessor = get(v, "preprocessor", "")
                
                var_stat = String.(split(k, "_"))
                statistic = length(var_stat) == 2 && var_stat[1] == variable ? var_stat[2] : nothing
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    # TODO: when would this happen??
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = MetaData(
                    variable, experiment, alias; 
                    timerange = timerange,
                    subdir = k,
                    variable_long = var_long,
                    statistic = statistic,
                    esmvaltoolPreproc = preprocessor
                )
                # if there are different recipes for CMIP5 and CMIP6, 'meta' might already 
                # have been added to 'meta_attribs'
                present = any(x -> isEqual(meta, x), meta_attribs)
                if !present
                    push!(meta_attribs, meta)
                end
            end
        end
    end
    if !isnothing(constraint)
        constrainMetaData!(meta_attribs, constraint)
    end
    return meta_attribs
end


function isEqual(m1::MetaData, m2::MetaData)
    return all(x -> getfield(m1, x) == getfield(m2, x), fieldnames(MetaData))
end


function metadataToDict(meta::MetaData; exclude::Vector{Symbol}=Vector{Symbol}())
    d = Dict{String, Any}()
    for f in filter(x -> !(x in exclude), fieldnames(MetaData))
        d[String(f)] = getfield(meta, f)
    end
    return d
end


function structToDict(data)
    return Dict(
        name => getfield(data, name) for name in fieldnames(typeof(data))
    )
end


function alignTimeseries!(data::Vector{<:YAXArray})
    if length(unique(map(x -> size(otherdims(x, :time)), data))) != 1
        throw(ArgumentError("Dimension sizes must be identical across data to align timeseries!"))
    end
    if !all(map(x -> hasdim(x, :time), data))
        throw(ArgumentError("All datasets must have time dimension to align timeseries!"))
    end
    
    year_min = minimum(map(x -> minimum(map(Dates.year, dims(x, :time))), data))
    year_max = maximum(map(x -> maximum(map(Dates.year, dims(x, :time))), data))
    nb_years = year_max - year_min + 1
    timerange = DateTime(year_min):Year(1):DateTime(year_max)
    
    none_time_dims = otherdims(data[1], :time)
    s = map(length, none_time_dims)
    T = Union{Missing, Float32}
    for (i, ds) in enumerate(data)
        # if ds allows missing values, undef is initialized with missing
        dat = Array{T}(undef, s..., nb_years)
        ds_extended = YAXArray((none_time_dims..., Dim{:time}(timerange)), dat, ds.properties)
        ds_extended[time = Where(x -> Dates.year(x) in map(Dates.year, dims(ds, :time)))] = ds
        data[i] = ds_extended    
    end
end


function overlapTimeseries!(data::Vector{AbstractArray})
    if !all(map(x -> hasdim(x, :time), data))
        throw(ArgumentError("All datasets must have time dimension to align timeseries!"))
    end
    start_years = map(x -> minimum(map(Dates.year, dims(x, :time))), data)
    end_years =  map(x -> maximum(map(Dates.year, dims(x, :time))), data)
    start_y = maximum(start_years)
    end_y = minimum(end_years)
    if start_y > end_y
        ts = unique(zip(start_years, end_years))
        throw(ArgumentError("The different timeseries to be merged do not overlap! Found timeseries $ts, latest start: $start_y  earlist end: $end_y"))
    end
    timerange = DateTime(start_y):Year(1):DateTime(end_y)
    timerange_years = map(Dates.year, timerange)
    for (i, ds) in enumerate(data)
        ds_years = map(x -> Dates.year(x), dims(ds, :time))
        indices = findall(x -> x in timerange_years, ds_years)
        ds_slice = ds[time = indices]
        data[i] = YAXArray((otherdims(ds, :time)..., Dim{:time}(timerange)), ds_slice.data, ds.properties)   
    end
end


"""
    metaDataFromYAML(ds::Dict)

Return metadata of data specified in `ds`.

# Arguments:
- `ds`: content of config yaml file specifying meta attributes and paths of data.
"""
function metaDataFromYAML(ds::Dict)
    fn_err(x) = throw(ArgumentError("$(x) must be provided for each dataset in config yaml file!"))
    experiment = get(() -> fn_err("exp"), ds, "exp")
    variables = get(() -> fn_err("variables"), ds, "variables")
    aliases = get(() -> fn_err("aliases"), ds, "aliases")
    
    subdirs = get(ds, "subdirs", nothing)
    statistics = get(ds, "statistics", nothing)

    meta_ds = Vector{MetaData}()
    if isnothing(subdirs)
        subdirs = isnothing(statistics) ? 
            fn_err("if subdirs not provided, statistics") : 
            combineAll(variables, statistics)
        for alias in aliases
            for subdir in subdirs
                var, stats = string.(split(subdir, "_"))
                push!(meta_ds, MetaData(var, experiment, alias; subdir=subdir, statistic=stats))     
            end
        end
    else
        for var in variables
            for alias in aliases
                for subdir in subdirs
                    push!(meta_ds, MetaData(var, experiment, alias; subdir=subdir))
                end
            end
        end
    end
    return meta_ds
end


function constrainSubdirs!(paths::Vector{String}, subdir_constraints::Vector{String})
    if !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), paths)
    end
    return nothing
end


"""
    resolvePathsFromMetaData( 
        meta::MetaData,
        base_path::String,
        dir_per_var::Bool;
        constraint::Union{Dict, Nothing} = nothing
    )

Return paths to data files for data specified in `meta`, possibly constraint by values in 
`constraint`. The paths were the data is stored is expected to follow the following 
structure (corresponding to the output from ESMValTool used for preprocessing the data):

`base_path` is the top-level directory. If `dir_per_var` is true, `base_path` is assumed to 
have a (or several) subdirectory for each climate variable with _VAR as part of the 
subdirectory's name (e.g. _tas, cmip5_tas, etc.). These subdirectories may be constraint by 
containing at least one of the values in `constraint["base_subdirs"]`. 

Let BASE refer to `base_path`, or respectively, to the subdirectories for the climate 
variables. Then the following structure is: BASE/preproc/meta.alias/meta.subdir. 
In ESMValTool, `meta.alias` corresponds to the (self-chosen) name under the section 
'diagnostics' and `meta.subdir` to the (self-chosen) name under the section 'variables'.

The returned paths are the paths to all files within this directory, possibly constraint by 
the filenames containig at least one string in `constraint["projects"]` and respectively 
at least one string in `constraint["models"]`. 
"""
function resolvePathsFromMetaData(
    meta::MetaData, 
    base_path::String, 
    dir_per_var::Bool;
    constraint::Union{Dict, Nothing} = nothing
)
    has_constraint = !isnothing(constraint) && !isempty(constraint)
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join = true)) # all subdirectories
        filter!(x -> occursin("_" * meta.variable, x), base_paths) # just subdirs for variable
        if isempty(base_paths)
            throw(ArgumentError("$(base_path) doesnt contain directories for variable $(meta.variable)!"))
        end
        if has_constraint
            bps = copy(base_paths)
            constrainSubdirs!(base_paths, get(constraint, "base_subdirs", Vector{String}()))
            if isempty(base_paths)
                throw(ArgumentError("$(bps) dont match with constraint base_subdirs: $(constraint["base_subdirs"])!"))
            end
        end
    else
        base_paths = [base_path]
    end
    paths_data_dirs_all = map(base_paths) do p
        joinpath(p, "preproc", meta.alias, meta.subdir)
    end
    paths_data_dirs = filter(isdir, paths_data_dirs_all)
    if isempty(paths_data_dirs)
        throw(ArgumentError("$paths_data_dirs_all arent't existing directories!"))
    end
    paths_to_files = map(p -> collectNCFilePaths(p), paths_data_dirs)
    return vcat(paths_to_files...)
end

"""
    parseFilename(filename::String, format::String)

Retrieve information from `filename` and return it as an object of type FilenameMeta.
"""
function parseFilename(filename::String, format::String)
    parts_format = split(format, "_")
    parts_fn = split(basename(splitext(filename)[1]) , "_")
    if length(parts_format) != length(parts_fn)
        throw(ArgumentError("File $filename doesn't align with format $format !"))
    end
    d = Dict{Symbol, String}()
    for k in fieldnames(FilenameMeta)
        indices = findall(x -> Symbol(lowercase(x)) == k, parts_format)
        if length(indices) > 1
            throw(ArgumentError("Invalid filename format: $format, $x appears more than once!"))
        elseif !isempty(indices)
            d[k] = parts_fn[indices[1]]
        end
    end
    return FilenameMeta(
        variable = get(d, :variable, ""),
        tableid = get(d, :tableid, ""),
        model =  get(d, :model, ""),
        experiment = get(d, :experiment, ""),
        variant = get(d, :variant, ""),
        fn = filename,
        grid = get(d, :grid, ""),
        timerange = get(d, :timerange, ""),
        mip = get(d, :mip, "")
    )
end

function parseFilename(filename::String, format::Symbol)
    if format == :esmvaltool
        n = length(split(filename, "_"))
        format = n == length(split(ESMVT_FORMAT_CMIP5, "_")) ? :esmvaltool_cmip5 : :esmvaltool_cmip6
    end
    err_msg = "Only filename formats $(keys(FILENAME_FORMATS)) and :esmvaltool are defined. Found: $format."
    fn_format = get(() -> throw(ArgumentError(err_msg)), FILENAME_FORMATS, format)
    return parseFilename(filename, fn_format)
end


function isRetained(fn_meta::FilenameMeta, constraint::Dict)
    fn_err(k::Symbol) = throw(ArgumentError("$k is not a valid metadata key for filenames. Allowed are: $(fieldnames(FilenameMeta))"))
    isOk(meta::Dict, key_meta::Symbol, constraint_vals::Vector{String}) = begin
        keep = true
        if !isempty(constraint_vals)
            value = get(() -> fn_err(key_meta), meta, key_meta)
            grid_val = get(() -> fn_err(:grid), meta, :grid)
            variant_val = get(() -> fn_err(:variant), meta, :variant)
            if isempty(value)
                keep = true
            else
                # Model constraint can refer to model or member (given as model#variant)
                if key_meta == :model
                    constraint_mms = map(x -> string.(split(x, MODEL_MEMBER_DELIM)), constraint_vals)
                    constraint_models = map(x -> x[1], constraint_mms)
                    constraint_members = map(constraint_mms) do x
                        length(x) == 2 ? string.(split(x[2], "_")) : [nothing]
                    end
                    constraint_variants = map(x -> x[1], constraint_members)
                    # for cmip6 models we added the grids to the member id (model#variant_grid)
                    constraint_grids = map(constraint_members) do x
                        length(x) == 2 ? x[2] : nothing
                    end
                    values_ok = map(constraint_models, constraint_variants, constraint_grids) do model, variant, grid
                        value == model && (absent(variant) || variant_val == variant) && (absent(grid) || grid_val == grid)
                    end
                    keep = any(values_ok)
                elseif key_meta == :fn
                    keep = any(map(x -> occursin(x, value), constraint_vals))
                else
                    keep = value in constraint_vals
                end
            end
        end
        if !keep
            @debug "File with meta info $fn_meta excluded due to constraint $key_meta ($constraint_vals)"
        end
        return keep
    end
    meta = structToDict(fn_meta)
    return isOk(meta, :variable, get(constraint, "variables", Vector{String}())) &&
        isOk(meta, :model, get(constraint, "models", Vector{String}())) &&
        isOk(meta, :variant, get(constraint, "variants", Vector{String}())) &&
        isOk(meta, :experiment, get(constraint, "experiments", Vector{String}())) &&
        isOk(meta, :mip, get(constraint, "mips", Vector{String}())) &&
        isOk(meta, :timerange, get(constraint, "timeranges", Vector{String}())) &&
        isOk(meta, :grid, get(constraint, "grids", Vector{String}())) &&
        isOk(meta, :tableid, get(constraint, "table_ids", Vector{String}())) &&
        isOk(meta, :fn, get(constraint, "filenames", Vector{String}()))
end


function collectNCFilePaths(path_data_dir::String)
    if !isdir(path_data_dir)
        throw(ArgumentError("$path_data_dir is not a directory!"))
    end
    paths_to_files = filter(x -> isfile(x) && endswith(x, ".nc"), readdir(path_data_dir; join=true))
    if isempty(paths_to_files)
        @warn "No .nc files found!"
    end
    return paths_to_files
end


"""
    constrainMetaData!(meta_attributes::Vector{MetaData}, constraint::Dict)

Subset entries in `meta_attributes` so that only those with properties specified 
in `constraint` remain.

# Arguments
- `constraint::Dict`: Mapping to vector specifiying the properties of which at least one 
must be present for an id to be retained.
"""
function constrainMetaData!(meta_attributes::Vector{MetaData}, constraint::Dict)
    timerange_constraints = get(constraint, "timeranges", Vector{String}())
    alias_constraints = get(constraint, "aliases", Vector{String}())
    timerangeOk(attrib::MetaData) = any(x -> attrib.timerange == x, timerange_constraints)
    aliasOk(attrib::MetaData) = any(x -> attrib.alias == x, alias_constraints)

    # timerange and alias don't have to match, it's sufficient if either timerange or alias match
    if !isempty(timerange_constraints) && !isempty(alias_constraints)
        filter!(x -> timerangeOk(x) || aliasOk(x), meta_attributes)
        if isempty(meta_attributes)
            msg = "Neither timeranges $(timerange_constraints) nor aliases $(alias_constraints) found in data!"
            throw(ArgumentError(msg))
        end
    elseif !isempty(timerange_constraints)
        filter!(timerangeOk, meta_attributes)
    elseif !isempty(alias_constraints)
        filter!(aliasOk, meta_attributes)
    end
    stats_constraints = get(constraint, "statistics", Vector{String}())
    vars_constraints = get(constraint, "variables", Vector{String}())
    stats_ok(attrib::MetaData) = isempty(stats_constraints) || any(x -> attrib.statistic == x, stats_constraints)
    vars_ok(attrib::MetaData) = isempty(vars_constraints) || any(x -> attrib.variable == x, vars_constraints)
    filter!(stats_ok, meta_attributes)
    filter!(vars_ok, meta_attributes)
    return nothing
end

"""
    avgObsDatasets(observations::YAXArray)

Take mean across dimension 'model'. Returned YAXArray does not have dimension 'model'.
"""
function avgObsDatasets(observations::YAXArray)
    if !hasdim(observations, :model)
        throw(ArgumentError("Obs data requires dim :model when averaging across datasets."))
    end
    if length(dims(observations, :model)) > 1
        return mean(observations, dims=:model)[model=At("combined")]
    else
        return observations[model=1]
    end
end


function distancesData(model_data::DataMap, obs_data::DataMap, config::Dict{String, Number})
    return distancesData(model_data, obs_data,  activeDiagnostics(config))
end

"""
    distancesData(model_data::DataMap, obs_data::DataMap, diagnostics::Vector{String})

Compute RMSEs between models and observations for `diagnostics`.

# Arguments:
- `diagnostics::Vector{String}`: keys for which values must be provided in `model_data` and 
`obs_data`.
"""
function distancesData(model_data::DataMap, obs_data::DataMap, diagnostics::Vector{String})
    return map(diagnostics) do d
        model = get(model_data, d, nothing)
        obs = get(obs_data, d, nothing)
        if isnothing(model) || isnothing(obs)
            throw(ArgumentError("Data missing for diagnostic $d."))
        end
        distancesData(model, obs)
    end
end

"""
    distancesData(models::YAXArray, observations::YAXArray)

Compute the distance as the area-weighted RMSE between model predictions and observations.

If `observations' has dimension 'model' (for several observational datasets), the average across all is taken.

# Arguments:
- `models::YAXArray`: dimensions must be 'lon', 'lat', 'model'/'member'.
- `observations::YAXArray`: dimensions must be 'lon', 'lat' and possibly 'model'.
- `metric::Symbol`: :rmse for Root Mean Squared Error (default) or :mse for Mean Squared Error
"""
function distancesData(models::YAXArray, observations::YAXArray; metric::Symbol=:rmse)
    throwErrorIfDimMissing(models, [:lon, :lat]; include = :all)
    obs = hasdim(observations, :model) ? avgObsDatasets(observations) : observations
    model_dim = modelDim(models)
    model_names = lookup(models, model_dim)
    if otherdims(models, model_dim) != dims(obs)
        msg = "Found: Model dims (other than model dim): $(DimensionalData.name.(otherdims(models, model_dim))) and observational dims: $(DimensionalData.name.(dims(obs)))"
        throw(ArgumentError("dimensions of model (other than model dim) and observation data must align!" * msg))
    end
    latitudes = collect(lookup(models, :lat))
    obs_data = collect(obs)
    model_data = collect(models)
    mask_missing =  (ismissing.(obs_data) .+ ismissing.(model_data)) .> 0 # observations or model is missing (or both)
    # mask_missing has same dimensions as model_data (due to broadcasting)
    idx_model = dimnum(models, model_dim)
    distances = similar(model_names, eltype(model_data))
    if !(metric in [:rmse, :mse])
        throw(ArgumentError("Optional argument :metric must be :rmse (default) or :mse! Found: $metric"))
    end
    fn = metric == :rmse ? areaWeightedRMSE : areaWeightedMSE
    @inbounds for i in eachindex(model_names)
        data_m = selectdim(model_data, idx_model, i)
        mask = selectdim(mask_missing, idx_model, i)
        aw_mat = areaWeightMatrix(latitudes, mask)
        distances[i] = fn(data_m, obs_data, aw_mat)
    end
    return YAXArray((Dim{model_dim}(model_names),), distances)
end



"""
    distancesData(
        model_data::AbstractArray, 
        obs_data::AbstractArray,
        model_names::AbstractArray{<:String},
        latitudes::AbstractArray{<:Real}; 
        metric::Symbol=:rmse
    )

Compute the distance as the area-weighted RMSE (default) or MSE between model predictions and observations.

# Arguments:
- `model_data::AbstractArray`: dimensions must be 'lon', 'lat' and 'model'/'member' in 3rd dimension.
- `obs_data::AbstractArray`: dimensions must be 'lon', 'lat'.
- `latitudes::AbstractArray{<:Real}`: values of latitudes in `model_data` and `obs_data`. 
- `metric::Symbol`: :rmse for Root Mean Squared Error (default) or :mse for Mean Squared Error,
- `idx_model::Int`: number of model/member-dimension in `model_data`.
"""
function distancesData(
    model_data::AbstractArray, 
    obs_data::AbstractArray, 
    latitudes::AbstractArray{<:Real}; 
    metric::Symbol=:rmse,
    idx_model::Int = 3
)
    s = size(model_data)
    if length(s) != length(size(obs_data))
        obs_data = insertSingletonDim(obs_data, idx_model)
    end
    mask_missing =  (ismissing.(obs_data) .+ ismissing.(model_data)) .> 0 # observations or model is missing (or both)
    obs_data = selectdim(obs_data, idx_model, 1)
    
    distances = Vector{eltype(model_data)}(undef, size(model_data)[idx_model])
    if !(metric in [:rmse, :mse])
        throw(ArgumentError("Optional argument :metric must be :rmse (default) or :mse!"))
    end
    fn = metric == :rmse ? areaWeightedRMSE : areaWeightedMSE
    @inbounds for i in 1:s[idx_model]
        data_m = selectdim(model_data, idx_model, i)
        mask = selectdim(mask_missing, idx_model, i)
        aw_mat = areaWeightMatrix(latitudes, mask)
        distances[i] = fn(data_m, obs_data, aw_mat)
    end
    return distances
end

"""

    distancesData(
        models::AbstractArray, observations::AbstractArray, latitudes::AbstractVector
    )

Compute the distance as the area-weighted RMSE between model predictions and observations.

Neither the observations nor the model data must contain missing values.

# Arguments:
- `models`: dimensions must be 'lon', 'lat', 'model'/'member' (in this order)
- `observations`: dimensions must be 'lon', 'lat' (in this order)
- `latitudes`: for computing area weight matrices
"""
function distancesData(
    models::AbstractArray, observations::AbstractArray, latitudes::AbstractVector
)
    if (size(models, 1) != size(observations, 1)) || (size(models, 2) != size(observations, 2))
        # throw(ArgumentError("Dimensions mismatch between model ($(size(models))) and data ($(size(data)))!"))
        throw(ArgumentError("Dimensions mismatch between model and data!"))
    end
    mask = areaWeightMatrix(latitudes, fill(false, size(models)))
    squared_diff = (models .- observations) .^ 2
    mse = sum(mask .* squared_diff, dims=(1,2))
    return sqrt.(vec(mse))
end




function distancesModels(model_data::DataMap, config::Dict{String, Number})
    return distancesModels(model_data, activeDiagnostics(config))
end

"""
    distancesModels(model_data::DataMap, diagnostics::Vector{String})

Compute the model-model distances for `model_data` at every diagnostic in `diagnostics` and 
return a vector of YAXArrays.
"""
function distancesModels(model_data::DataMap, diagnostics::Vector{String})
    return map(diagnostics) do d
        model = get(model_data, d, nothing)
        if isnothing(model)
            throw(ArgumentError("Data missing for diagnostic $d."))
        end
        distancesModels(model)
    end
end

"""
    distancesModels(data::YAXArray)

Compute the area weighted RMSE between model predictions for each pair of models.

# Arguments:
- `data::YAXArray`: first two dimensions must be 'lon', 'lat', third 'model' or 'member'.
"""
function distancesModels(data::YAXArray)
    dim_symbol = modelDim(data)
    models = collect(lookup(data, dim_symbol))
    idx_model = dimnum(data, dim_symbol)
    nb_models = length(models)
    
    arr = Array(data)
    # only take values where none (!) of the models have missing values!! (Not just the two that are compared to one another)
    mask_missing = dropdims(any(ismissing, arr, dims=idx_model), dims = idx_model)
    aw_mat = areaWeightMatrix(collect(lookup(data, :lat)), mask_missing)
    matrixS = zeros(nb_models, nb_models)
    @inbounds for i in range(1, length(models)-1)
        @inbounds for j in range(i+1, length(models))
            m1 = selectdim(arr, idx_model, i)
            m2 = selectdim(arr, idx_model, j)
            matrixS[i, j] = areaWeightedRMSE(m1, m2, aw_mat)
        end
    end
    symDistMatrix = matrixS .+ matrixS'
    new_dim1 = dim_symbol == :member ? :member1 : model1
    new_dim2 = dim_symbol == :member ? :member2 : model2
    return YAXArray(
        (Dim{new_dim1}(models), Dim{new_dim2}(models)), symDistMatrix, deepcopy(data.properties)
    )
end


"""
    generalizedDistances(
        distances_all::YAXArray, weights::YAXArray; norm_avg_members::Bool = true
    )

For every diagnostic in `distances_all`, compute the weighted sum of all diagnostics.

# Arguments:
- `distances_all::YAXArray`: must have dimension :diagnostic.
- `weights::YAXArray`: weights assigned to each diagnostic.
- `norm_avg_members::Bool`: if true (default), average distances per model are computed 
BEFORE computing the median that is used to normalize the distances (seperately for each 
diagnostic); if false, median is computed across all individual members of each diagnostic.
"""
function generalizedDistances(distances_all::YAXArray, weights::YAXArray; norm_avg_members::Bool = true)
    throwErrorIfDimMissing(distances_all, :diagnostic)
    model_dims = modelDims(distances_all)
    norm = dropdims(median(distances_all, dims=model_dims), dims=model_dims)

    has_model_dim = :model in model_dims
    if !has_model_dim
        distances_all = :member1 in  model_dims ?
            summarizeMembersMatrix(distances_all, false) :
            summarizeMembersVector(distances_all)
    end
    dimensions = modelDims(distances_all)
    norm = norm_avg_members ? dropdims(median(distances_all, dims=dimensions), dims=dimensions) : norm
    normalized_distances = @d distances_all ./ norm
    
    weighted_dists = @d normalized_distances .* weights
    return dropdims(sum(weighted_dists, dims = :diagnostic), dims=:diagnostic)
end


"""
    generalizedDistances(
        distances_all::Vector{<:YAXArray}, diagnostics::Vector{String}, weights::YAXArray
    )

# Arguments:
- `distances_all::Vector{<:YAXArray}`: each entry refers to the distances of one diagnostic, 
each must have dimension 'model', 'member' or 'member1' and 'member2' (must be identical for 
every entry!).
- `diagnostics:: Vector{String}`: names of diagnostics.
- `weights::YAXArray`: weights assigned to each diagnostic.
- `norm_avg_members::Bool`: if true (default), average distances per model are computed 
BEFORE computing the median that is used to normalize the distances (seperately for each 
diagnostic); if false, median is computed across all individual members of each diagnostic.
"""
function generalizedDistances(
    distances_all::Vector{<:YAXArray}, diagnostics::Vector{String}, weights::YAXArray;
    norm_avg_members::Bool = true
)
    if length(distances_all) != length(diagnostics)
        throw(ArgumentError("distances must be computed for every diagnostic, found $(length(distances_all)) distances, but $(length(diagnostics)) diagnostics."))
    end
    map(x -> throwErrorIfDimMissing(x, [:model, :member, :member1]; include=:any), distances_all)
    normalizations = median.(distances_all)
    # average model members
    model_dims = modelDims.(distances_all)
    has_model_dim = all(t -> :model in t, model_dims)
    if !has_model_dim
        distances_all = all(t -> :member1 in t, model_dims) ?
            summarizeMembersMatrix.(distances_all, false) :
            summarizeMembersVector.(distances_all)
    end
    normalizations = norm_avg_members ? median.(distances_all) : normalizations
    normalized_distances = map(distances_all, normalizations, diagnostics) do dists, norm, diagnostic
        d = Array(dists) ./ norm
        dimensions = dims(dists)
        YAXArray(
            (dimensions..., Dim{:diagnostic}([diagnostic])), 
            reshape(d, (size(dimensions)..., 1)), 
            deepcopy(dists.properties)
        )
    end
    distances = cat(normalized_distances...; dims=Dim{:diagnostic}(diagnostics))
    weighted_dists = @d distances .* weights
    return mapslices(sum, weighted_dists; dims = (:diagnostic,))
end




function fixModelNamesMetadata(names::Vector{String})
    model_names = copy(names)
    for m in keys(MODEL_NAME_FIXES)
        indices = findall(model_names .== m)
        if !isempty(indices)
            model_names[indices] .= MODEL_NAME_FIXES[m]
        end
    end
    return model_names
end


function fixModelNameMetadata(name::String)
    return get(MODEL_NAME_FIXES, name, name)
end


function getMask(orog_data::YAXArray; mask_out_land::Bool = true)
    ocean_mask = orog_data .== 0
    indices_missing = findall(x -> ismissing(x), ocean_mask)
    # load data into memory with Array() for modification
    ocean_mask_mat = Array(ocean_mask)
    ocean_mask_mat[indices_missing] .= false
    meta = deepcopy(orog_data.properties)
    meta["_ref_id_mask"] = orog_data.properties["id"]
    mask_arr = YAXArray(dims(ocean_mask), Bool.(ocean_mask_mat), meta)
    return mask_out_land ? mask_arr : mask_arr .== false
end


"""
    sharedLevelMembers(data::DataMap)

Return vector of model members that are shared across all entries in `data`.
"""
function sharedLevelMembers(data::DataMap)
    if !(all(((id, dat),) -> hasdim(dat, :member), data))
        throw(ArgumentError("All datasets must have dimension :member!"))
    end
    members = map(x -> dims(x, :member), collect(values(data)))
    return reduce(intersect, members)
end


"""
    sharedLevelModels(data::DataMap)

Return vector of models (on level of models, not members) that are shared across all entries in `data`.
"""
function sharedLevelModels(data::DataMap)
    all_models = Vector(undef, length(data))
    for (i, id) in enumerate(keys(data))
        ds = data[id]
        if hasdim(ds, :model)
            all_models[i] = dims(ds, :model)
        elseif hasdim(ds, :member)
            all_models[i] =
                map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), dims(ds, :member))
        else
            throw(ArgumentError("Data for $id must have dimension :member or :model!"))
        end
    end
    return reduce(intersect, all_models)
end


"""
    loadModelsFromCSV(
        path::String, col_models::String; col_variants::Union{String,Nothing}=nothing
    )

Return a vector with models retrieved from csv file at `path`. If col_variants is provided, 
returned models are on level of model members (MODEL#variant, e.g. AWI-ESM#r1i1p1f1).
"""
function loadModelsFromCSV(
    path::String, col_models::String; col_variants::Union{String, Nothing}=nothing
)
    models = DataFrame(CSV.File(path))
    dropmissing!(models, col_models)
    ids = models[!, col_models]
    if !isnothing(col_variants)
        variants = map(s -> String.(strip.(split(s, ";"))), models[!, col_variants])
        for (i, m) in enumerate(models[!, col_models])
            variants[i] = map(v -> join([m, v], MODEL_MEMBER_DELIM), variants[i])
        end
        ids = unique(vcat(variants...))
    end
    return String.(ids)
end


"""
    areaWeightedMSE(m1::AbstractArray, m2::AbstractArray, aw_mat::AbstractArray)

Compute the area weighted (approximated by cosine of latitudes in radians) mean squared 
error between `m1` and `m2`.

# Arguments:
- `m1`: must have dimensions 'lon', 'lat' as first dimensions.
- `m2`: must have dimensions 'lon', 'lat' as first dimensions.
- `aw_mat`: matrix with area weights, of same size as `m1` and `m2`, that will be normalized.
"""
function areaWeightedMSE(m1::AbstractArray, m2::AbstractArray, aw_mat::AbstractArray)
    if size(m1) != size(m2) || size(m1) != size(aw_mat)
        throw(ArgumentError("All input arrays must have same size to compute areaweighted rmse! Found: $(size.([m1, m2, aw_mat]))."))
    end
    squared_diff = (m1 .- m2) .^ 2
    return sum(skipmissing(aw_mat .* squared_diff))./ sum(aw_mat)
end



"""
    areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, aw_mat::AbstractArray)

Compute the area weighted (approximated by cosine of latitudes in radians) root mean squared 
error between `m1` and `m2`. 

# Arguments:
- `m1`: must have dimensions 'lon', 'lat' as first dimensions.
- `m2`: must have dimensions 'lon', 'lat' as first dimensions.
- `aw_mat`: matrix with area weights, of same size as `m1` and `m2`, that will be normalized.
"""
function areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, aw_mat::AbstractArray)
    return sqrt(areaWeightedMSE(m1, m2, aw_mat))
end

"""
    areaWeightedSquaredErr()

# Arguments:
- `model_data`: must have same dimensions as `obs_data`
- `obs_data`: must have same dimensions as `model_data`
- `latitudes`: 
"""
function areaWeightedSquaredErr(
    model_data::AbstractArray, 
    obs_data::AbstractArray, 
    latitudes::AbstractArray{<:Real}; 
)
    mask =  (ismissing.(obs_data) .+ ismissing.(model_data)) .> 0 # observations or model is missing (or both)
    aw_mat = areaWeightMatrix(latitudes, mask)
    squared_diff = (model_data .- obs_data) .^ 2
    return skipmissing(aw_mat .* squared_diff) ./ sum(aw_mat)
end



"""
    uncertaintyRanges(
        data::AbstractArray; w::Union{DimArray, Nothing}=nothing, quantiles=[0.167, 0.833]
    )

Compute weighted `quantiles` of timeseries `data`.

# Arguments:
- `data::AbstractArray`: must have dimensions 'time', 'model'.
- `w::Union{DimArray, Nothing}=nothing`: must have dimensions 'model' and 'weight', each 
'weight'-vector must sum up to 1.
- `quantiles=[0.167, 0.833]`: vector with two entries between 0 and 1 representing  the 
lower and upper bound in this order.
"""
function uncertaintyRanges(
    data::YAXArray; 
    weights::Union{YAXArray, Nothing} = nothing, 
    q_lower::Number = 0.167, q_upper::Number = 0.833
)
    model_dim = modelDim(data)
    meta = deepcopy(data.properties)
    meta["_confidence_quantiles"] = string.([q_lower, q_upper])
    # TODO: check missing values!!
    #     if sum(ismissing.(slice)) >= length(slice) - 1 # at least two values must be given
    #         qs = [missing, missing]
    #     end
    if isnothing(weights) || !hasdim(weights, :weight)
        if !isnothing(weights)
            weights, data = alignWeightsAndData(data, weights)
        end
        qs = mapslices(
            slice -> map(q -> quantile(coalesce.(slice, NaN), q; w=weights), [q_lower, q_upper]),
            data,
            dims=(model_dim,)
        )
        qs = setDim(qs, :OutAxis1, :confidence, ["lower", "upper"])
    else
        quantiles = []
        for weight_id in lookup(weights, :weight)
            w, df = alignWeightsAndData(data, weights[weight = At(weight_id)])
            qs = mapslices(
                slice -> map(q -> quantile(coalesce.(slice, NaN), q; w=w), [q_lower, q_upper]),
                df,
                dims=(model_dim,)
            )
            qs = setDim(qs, :OutAxis1, :confidence, ["lower", "upper"])
            push!(quantiles, qs)
        end
        qs = concatenatecubes(quantiles, dims(weights, :weight))
    end
    return YAXArray(dims(qs), qs.data, meta)
end


"""
    quantile(samples::Vector{<:Number}, p::Number)

Compute the p-th quantile.
"""
function quantile(
    samples::AbstractArray{<:Number}, 
    p::Number;
    w::Union{Nothing, AbstractArray{<:Number}} = nothing,
    interpolate::Bool = true
)
    if p < 0 || p > 1
        throw(ArgumentError("It must hold: 0<p<=1! Found: $p"))
    end
    indices = sortperm(samples)
    samples = samples[indices]
    n = length(samples)
    if isnothing(w)
        probs_cum = cumsum(repeat([1/n], n))
    else
        probs_cum = cumsum(w[indices])
    end
    idx = findfirst(x -> x >= p, probs_cum)
    
    h = (n - 1) * p + 1
    interpolate = isinteger(h) ? false : interpolate
    
    if interpolate
        if isnothing(w)
            # use linear interpolation default definition used in numpy
            m = 1 - p
            j = Int(floor(n * p + m))
            gamma = n * p + m - j
            result = (1-gamma) * samples[j] + gamma * samples[j+1]
        else
            # for interpolation + weighted, use empirical cdf:
            q1 = probs_cum[idx]
            if q1 == p
                result = samples[idx]
            else
                if idx==1
                    result = samples[1]
                else
                    # standard linear interpolation
                    x0, x1 = samples[idx - 1], samples[idx]
                    q0, q1 = probs_cum[idx - 1], probs_cum[idx]
                    result =  x0 + ((p - q0)/(q1 - q0)) * (x1 - x0)
                end
            end
        end
    else
        result = samples[idx]
        # without interpolation: this is the same
        # fn = StatsBase.ecdf(samples; weights=w)
        # probs = fn.(samples)
        # idx = findfirst(x -> x >=p, probs)
    end
    return result
end


"""
    metaDataChecksCMIP(meta::Dict, path::String)

Check model names as retrieved from the metadata for potential inconsistencies wrt filename.
"""
function metaDataChecksCMIP(meta::Dict{String, T}, path::String) where T <: Any
    name = meta[getCMIPModelsKey(Dict(meta))]
    if !occursin(name, basename(path)) && !(name in keys(MODEL_NAME_FIXES))
        @warn "model name as read from metadata of stored .nc file ($name) and used as dimension name is not identical to name appearing in filename $(basename(path))."
    end 
    return nothing
end


function maskFileConstraints(
    paths::Vector{String}, fn_format::Union{Symbol, String}, constraint::Dict
)
    filenames = first.(splitext.(basename.(paths)))
    filenames_meta = parseFilename.(filenames, fn_format)
    return isRetained.(filenames_meta, fill(constraint, length(filenames_meta)))
end


"""
    modelDim(data::YAXArray)

Return the model dimension of `data` which must either be :model or :member.

If none is present, throw ArgumentError.
"""
function modelDim(data::YAXArray)
    err_msg = "Data must have either dimension :member or :model, found: $(dims(data))."
    return hasdim(data, :model) ? :model : 
        hasdim(data, :member) ? :member : throw(ArgumentError(err_msg))
end

"""
    modelDims(data::YAXArray)

Return vector of dimensions of `data` that contain either 'model' or 'member'.
"""
function modelDims(data::YAXArray)
    dim_names = string.(dimNames(data))
    model_dims = filter(d -> occursin("model", d), dim_names)
    if isempty(model_dims)
        # check if it contains dimension with name 'member'
        model_dims = filter(d -> occursin("member", d), dim_names)
    end
    if isempty(model_dims)
        throw(ArgumentError("Data must contain a dimension with 'member' or 'model' in its name! found: $dim_names"))
    end
    return Symbol.(model_dims)
end



function apply!(
    dm::DataMap,
    fn::Function,
    args...; 
    ids::AbstractVector{T} = Vector{String}(),
    ids_new::AbstractVector{T} = Vector{String}(),
    kwargs...
) where T <: Union{String, Symbol}

    ids = isempty(ids) ? collect(keys(dm)) : ids
    ids_new = isempty(ids_new) ? ids : ids_new
    for (id, id_new) in zip(ids, ids_new)
        dm[id_new] = fn(dm[id], args...; kwargs...)
    end
    return nothing
end

"""

    apply(
        dm::DataMap,
        fn::Function,
        args...; 
        ids::AbstractVector{T} = Vector{String}(),
        ids_new::AbstractVector{T} = Vector{String}(),
        kwargs...
    )

Apply `fn` with positional arguments `args` and keyword arguments 
`kwargs` and return result as a new DataMap where keys are `ids_new` if provided, otherwise
same keys as in `dm` are used.


# Arguments:
- `dm::DataMap`: data.
- `fn::Function`: function to be applied.
- `args...`: positional arguments for `fn`.
- `ids::AbstractVector{T}=Vector{String}()`: keys for data on which `fn` is applied; if empty all keys of `dm` are used.
- `ids_new::AbstractVector{T}=Vector{String}()`: keys used in new datamap; if empty `ids` are used instead.
- `kwargs...`: keyword arguments for `fn`.
"""
function apply(
    dm::DataMap,
    fn::Function,
    args...; 
    ids::AbstractVector{T} = Vector{String}(),
    ids_new::AbstractVector{T} = Vector{String}(),
    kwargs...
) where T <: Union{String, Symbol}

    dm_new = DataMap()
    ids = isempty(ids) ? collect(keys(dm)) : ids
    ids_new = isempty(ids_new) ? ids : ids_new
    for (id, id_new) in zip(ids, ids_new)
        dm_new[id_new] = fn(dm[id], args...; kwargs...)
    end
    return dm_new
end


"""
    indicesTimeseries(times::Vector, constraint_ts::Dict)

Get indices of `times` that are exactly aligned with `start_y` and `end_y` given in `constraint_ts`.
If no values are given, return indices for entire vector `times`.

"""
function indicesTimeseries(times::Vector{DateTime}, constraint_ts::Dict)
    start_y = get(constraint_ts, "start_y", -Inf)
    end_y = get(constraint_ts, "end_y", Inf)
    indices_time = findall(t -> Dates.year(t) >= start_y && Dates.year(t) <= end_y, times)
    times = times[indices_time]
    # only use data that's exactly from start to end!
    start_wrong = start_y != -Inf && (!isempty(times) && start_y != minimum(map(Dates.year, times)))
    end_wrong   = end_y != Inf && (!isempty(times) && end_y != maximum(map(Dates.year, times)))
    return (isempty(times) || start_wrong || end_wrong) ? [] : indices_time
end


function alignWeightsAndData(data::YAXArray, weights::YAXArray)
    level_data = Data.modelDim(data)
    level_weights = Data.modelDim(weights)
    if level_data != level_weights
        if level_data == :member
            data = summarizeMembers(data)
        else
            weights = summarizeMembers(weights; fn=sum)
        end
    end
    models_data = collect(dims(data, level_data))
    models_weights = collect(dims(weights, level_data))
    data_no_weights = [model for model in models_data if !(model in models_weights)]
    if !isempty(data_no_weights)
        msg = "No weights were computed for follwoing models, thus not considered in the weighted average:"
        @warn msg data_no_weights
        # Only include data for which there are weights
        indices = findall(x -> !(x in data_no_weights), models_data)
        data = indexModel(data, (level_data,), indices)
    end
    # weights for which we don't have data
    weights_no_data = filter(x -> !(x in models_data), models_weights)
    if !isempty(weights_no_data)
        @warn "Weights were renormalized since data of models missing for which weights have been computed: $weights_no_data"
        # renormalize weights
        indices = findall(m -> m in dims(data, level_data), models_weights)
        weights = indexModel(weights, (level_data,), indices)
        weights = weights ./ sum(weights)
    end
    return (weights, data)
end


"""
    subsetDataMap(data::DataMap, ids::Vector{String})

Return new DataMap with data from `data` at `ids`.
"""
function subsetDataMap(data::DataMap, ids::Vector{String})
    arrs = Vector{YAXArray}(undef, length(ids))
    for (i, id) in enumerate(ids)
        df = data[id]
        arrs[i] = YAXArray(dims(df), df.data, deepcopy(df.properties))
    end
    return defineDataMap(arrs, ids)
end


"""
    reduceNbMembers(dat_members::YAXArray; max_n = 5)

For every model in `dat_members` that has more than `max_int` members, randomly sample 
`max_int` members. Return YAXArray with reduced number of members per model.
"""
function reduceNbMembers(dat_members::YAXArray; max_n::Int = 5)
    nb_members_dict = nbModelMembers(dat_members)
    members = collect(dat_members.member)
    df_large = filter(((k,v),) -> v > max_n, nb_members_dict)
    indices_members = []
    for model in collect(keys(df_large))
        indices_model = findall(x -> startswith(x, model * MODEL_MEMBER_DELIM), members);
        include_model = StatsBase.sample(indices_model, max_n, replace=false)
        push!(indices_members, include_model)
    end
    df_ok = filter(((k,v),) -> v <= max_n, nb_members_dict)
    for model in collect(keys(df_ok))
        indices_model = findall(x -> startswith(x, model * MODEL_MEMBER_DELIM), members);
        push!(indices_members, indices_model)
    end
    indices_members = sort(collect(Iterators.flatten(indices_members)))
    return dat_members[member = indices_members]
end


"""
# Arguments:
- `data::AbstractArray`: lon x lat x model
- `n_rep::Int`: number of repetitions for selected model
- `index::Int`: index of model to be repeated (default: 2nd)
"""
function repeatModel(data::AbstractArray, n_rep::Int; index::Int=2)
    if n_rep == 1
        return copy(data)
    end
    models = Array(data.model)
    n_models = length(models)
    m_rep = models[index]
    models_rep = [map(i -> m_rep * "#" * string(i), 1:n_rep)...]
    models_new = vec(hcat([models[1:index-1]..., models_rep..., models[index+1:end]...]))
    s = size(data)
    data_rep = YAXArray(
        (dims(data, :lon), dims(data, :lat), Dim{:model}(models_new)),
        rand(s[1], s[2], length(models) + n_rep - 1)
    )
    if index > 1
        data_rep[:,:,1:index-1] .= Array(data[model=1:index-1])
    end
    for i in 1:n_rep
        data_rep[:, :, index+i-1] .= Array(data[model=index])
    end
    if index < n_models
        data_rep[:, :, index + n_rep : n_models + n_rep - 1] .= Array(data[model = index + 1 : n_models])
    end
    return data_rep
end

function initYAX(dimensions)
    return YAXArray(dimensions, fill(missing, size(dimensions)...))
end

function emptyYAX(arr::YAXArray)
    return all(ismissing.(arr))
end

# TODO
# function warnIfModelConstraintNotFulfilled(
#     constraints::Vector{<:Dict{<:Any, <:Any}}, 
#     loaded_data::DataMap, 
#     ids::Vector{String}
# )
#     model_constraints = map(x -> get(x, "models", Vector{String}()), constraints)
#     if length(unique(model_constraints)) == 1 && !isempty(model_constraints[1]) 
#         model_constraints = model_constraints[1:1]
#     end
#     pattern = Regex("^([^" * MODEL_MEMBER_DELIM * "]+)" * MODEL_MEMBER_DELIM * "([^_]+)_?(.*)")
#     for (i, requested_models) in enumerate(model_constraints)
#         ds = loaded_data[ids[i]]
#         found_members_with_grid = lookup(ds, :member)
#         found_members = Vector{String}(undef, length(found_members_with_grid))
#         found_models = Vector{String}(undef, length(found_members_with_grid))
#         for (i, member) in enumerate(found_members_with_grid)
#             m = match(pattern, member)
#             found_models[i] = m.captures[1]
#             found_members[i] = m.captures[1] * MODEL_MEMBER_DELIM * m.captures[2]            
#         end
#         all_not_found = filter(x -> !(x in found_models) && !(x in found_members) && !(x in found_members_with_grid), requested_models)
#         models_not_found = modelsFromMemberIDs(all_not_found; uniq=true)
#         if !isempty(models_not_found)
#             @warn "For the following requested models no data was found: " models_not_found
#         end
#     end
#     return nothing
# end