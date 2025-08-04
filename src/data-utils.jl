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
    data = deepcopy(data)
    dim_names = dimNames(data)
    index_vec = map(dim_names) do name 
        name in model_dims ? indices : Colon()
    end
    data = data[index_vec...]
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
For timeseries data, the time dimension may cover different ranges. In that case the 
maximal timeseries is used and filled with NaN for missing values.  

The combined YAXArray has the additional dimension `new_dim` (default: :model) with `names`
as values. If `names` is not provided, default values 'model1', 'model2, etc. are used 
instead.

If `sorted` is true, model dimension of returned data is sorted alphabetically and the 
vector entries in the metadata dictionary of the returned array are also sorted accordingly.
"""
function combineModelsFromMultipleFiles(
    data::AbstractVector{<:AbstractArray}; 
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
        if !all(map(x -> hasdim(x, :time), data))
            msg = "Data does not have the same size across all models: $(data_sizes)"
            throw(ArgumentError(msg))
        else
            # if difference only in time, use maximal possible timeseries and add NaNs
            #alignTimeseries!(data)
            overlapTimeseries!(data)
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
    dimData = YAXArray(dimData.axes, dimData.data, meta_dict)
    if length(model_names) != length(unique(model_names))
        duplicates = unique([m for m in model_names if sum(model_names .== m) > 1])
        @warn "Some datasets appear more than once" duplicates
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


function modelsFromMemberIDs(members::AbstractVector{<:String}; uniq::Bool=false)
    models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), members)
    return uniq ? unique(models) : models
end

function modelsFromMemberIDs(data::YAXArray; uniq::Bool=false)
    throwErrorIfDimMissing(data, :member)
    return modelsFromMemberIDs(string.(dims(data, :member)); uniq)
end


"""
    approxAreaWeights(latitudes::Vector{<:Number})

Create a YAXArray with the cosine of `latitudes` which approximates the cell 
area on the respective latitude.
"""
function approxAreaWeights(latitudes::AbstractVector{<:Number})
    area_weights = cos.(deg2rad.(latitudes))
    return YAXArray((Dim{:lat}(latitudes),), area_weights)
end

"""
    makeAreaWeightMatrix

# Arguments:
- mask: first two dimensions must refer to lon, lat
"""
function makeAreaWeightMatrix(
    longitudes::AbstractVector{<:Number},
    latitudes::AbstractVector{<:Number};
    mask::Union{YAXArray, Nothing} = nothing
)
    area_weights = approxAreaWeights(latitudes)
    area_weighted_mat = ones(length(longitudes)) * area_weights.data' # matrix multipication
    if !isnothing(mask)
        s_mask = size(mask)
        if length(s_mask) > 2
            if length(s_mask) > 3 
                throw(ArgumentError("area weight matrix only implemented for lonxlat and lonxlatxmodel input"))
            end
            n_models = s_mask[3]
            area_weighted_mat = cat(fill(area_weighted_mat, n_models)...; dims=3)
        end
        area_weighted_mat[mask] .= 0
    end
    area_weighted_mat = area_weighted_mat ./ sum(area_weighted_mat, dims=(1,2))
    dimensions = isnothing(mask) ? (Dim{:lon}(longitudes), Dim{:lat}(latitudes)) : dims(mask)
    return YAXArray(dimensions, area_weighted_mat)
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
            if !isempty(dups) && !isempty(warn_msg)
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
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model" # Test this (why?)
        data = setDim(data, dim, new_dim_name, models)
    end
    return data
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
    if !all(map(x -> hasdim(x, :time), data))
        throw(ArgumentError("All datasets must have time dimension to align timeseries!"))
    end
    year_min = minimum(map(x -> minimum(map(Dates.year, dims(x, :time))), data))
    year_max = maximum(map(x -> maximum(map(Dates.year, dims(x, :time))), data))
    nb_years = year_max - year_min + 1

    timerange = DateTime(year_min):Year(1):DateTime(year_max)
    for (i, ds) in enumerate(data)
        none_time_dims = otherdims(ds, :time)
        s = map(length, none_time_dims)
        # if ds allows missing values, undef is initialized with missing
        dat = Array{eltype(ds)}(undef, s..., nb_years)
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
    metaDataFromYAML(content::Dict)

Return metadata of data specified in `content` possibly constrained by values in `arg_constraint`.

For constraints that are specified in `content` as well as in the `arg_constraint` argument, 
the values of the latter have precedence over the former. The constraints given
in the argument `arg_constraint` are applied to EVERY dataset specified in the config 
file.

# Arguments:
- `content`: content of config yaml file specifying meta attributes and paths of data
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

Retrieve information from `filename` and returns it as an object of type FilenameMeta.
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
        isOk(meta, :fn, get(constraint, "filename", Vector{String}()))
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

If several observational datasets are present, the average across all is taken.
"""
function distancesData(models::YAXArray, observations::YAXArray)
    obs = hasdim(observations, :model) ? avgObsDatasets(observations) : observations
    # copy of data needed?
    maskMissing = @d (ismissing.(obs) .+ ismissing.(models)) .> 0 # observations or model is missing (or both)
    return areaWeightedRMSE(models, obs,  maskMissing)
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
    # data = deepcopy(data) copy of data needed?
    # only take values where none (!) of the models has infinite values!! (Not just the two that are compared to one another)
    dim_symbol = modelDim(data)
    nb_models = length(dims(data, dim_symbol))
    mask_missing = dropdims(any(ismissing, data, dims = dim_symbol), dims = dim_symbol)

    matrixS = zeros(nb_models, nb_models)
    for (i, model_i) in enumerate(eachslice(data[:, :, 1:(end-1)]; dims = dim_symbol))
        for (j, model_j) in enumerate(eachslice(data[:, :, (i+1):end]; dims = dim_symbol))
            idx = j + i
            matrixS[i, idx] = areaWeightedRMSE(model_i, model_j, mask_missing)
        end
    end
    symDistMatrix = matrixS .+ matrixS'
    dim = Array(dims(data, dim_symbol))
    new_dim1 = dim_symbol == :member ? :member1 : model1
    new_dim2 = dim_symbol == :member ? :member2 : model2
    return YAXArray((Dim{new_dim1}(dim), Dim{new_dim2}(dim)), symDistMatrix)
end


"""
    generalizedDistances(distances_all::YAXArray, weights::YAXArray)

For every diagnostic in `distances_all`, compute the weighted sum of all diagnostics.

# Arguments:
- `distances_all::YAXArray`: must have dimension :diagnostic.
"""
function generalizedDistances(distances_all::YAXArray, weights::YAXArray)
    throwErrorIfDimMissing(distances_all, :diagnostic)
    model_dims = modelDims(distances_all)
    has_model_dim = :model in model_dims
    if !has_model_dim
        distances_all = :member1 in  model_dims ?
            summarizeMembersMatrix(distances_all, false) :
            summarizeMembersVector(distances_all)
    end
    
    dimensions = modelDims(distances_all)
    norm = dropdims(median(distances_all, dims=dimensions), dims=dimensions)
    normalized_distances = @d distances_all ./ norm

    # normalized_distances = model_vs_data ? 
    #     summarizeMembersVector(normalized_distances) :
    #     summarizeMembersMatrix(normalized_distances, false)
    
    weighted_dists = @d normalized_distances .* weights
    return dropdims(sum(weighted_dists, dims = :diagnostic), dims=:diagnostic)
end


"""
    generalizedDistances(distances_all::Vector{YAXArray}, weights::YAXArray)

# Arguments:
- `distances_all::Vector{YAXArray}`: each entry refers to the distances of one diagnostic, 
each must have dimension 'model', 'member' or 'member1' and 'member2' (must be identical for 
every entry!).
- `weights::YAXArray`
"""
function generalizedDistances(
    distances_all::Vector{<:YAXArray}, diagnostics::Vector{String}, weights::YAXArray
)
    if length(distances_all) != length(diagnostics)
        throw(ArgumentError("distances must be computed for every diagnostic, found $(length(distances_all)) distances, but $(length(diagnostics)) diagnostics."))
    end
    map(x -> throwErrorIfDimMissing(x, [:model, :member, :member1]; include=:any), distances_all)

    # average model members
    model_dims = modelDims.(distances_all)
    has_model_dim = all(t -> :model in t, model_dims)
    if !has_model_dim
        distances_all = all(t -> :member1 in t, model_dims) ?
            summarizeMembersMatrix.(distances_all, false) :
            summarizeMembersVector.(distances_all)
    end
    # normalize
    normalizations = median.(distances_all)
    normalized_distances = map(distances_all, normalizations, diagnostics) do dists, norm, diagnostic
        d = @d dists ./ norm
        YAXArray(
            (dims(d)..., Dim{:diagnostic}([diagnostic])), 
            reshape(d.data, (size(dims(d))..., 1)), 
            d.properties
        )
    end
    # average; distances_all either has (in every entry) :member or :member1 and :member2
    # normalized_distances = all(hasdim.(distances_all, :member1)) ?
    #     summarizeMembersMatrix.(normalized_distances, false) :
    #     summarizeMembersVector.(normalized_distances)
    
    # for matrix, dimensions are symmetric, just take first, st. it also works for vector
    dimension = filter(x -> x != :diagnostic, dimNames(normalized_distances[1]))[1]
    models = map(x -> dims(x, dimension) , normalized_distances)
    shared_models = intersect(models...)
    
    distances = map(dists -> subsetModelData(dists, shared_models), normalized_distances)
    distances = cat(distances..., dims = Dim{:diagnostic}(diagnostics))
    weighted_dists = @d distances .* weights
    return dropdims(sum(weighted_dists, dims = :diagnostic), dims=:diagnostic)
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
        variants = map(s -> strip.(split(s, ";")), models[!, col_variants])
        for (i, m) in enumerate(models[!, col_models])
            variants[i] = map(v -> join([m, v], MODEL_MEMBER_DELIM), variants[i])
        end
        ids = vcat(variants...)
    end
    return String.(ids)
end


"""
    areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)

Compute the area weighted (approximated by cosine of latitudes in radians) root mean squared 
error between `m1` and `m2`. 

# Arguments:
- `m1`: must have dimensions 'lon', 'lat' as first dimensions.
- `m2`: must have dimensions 'lon', 'lat' as first dimensions.
- `mask`: indicates missing values (1 refers to missing, 0 to non-missing).
"""
function areaWeightedRMSE(m1::AbstractArray, m2::AbstractArray, mask::AbstractArray)
    if dims(m1, :lon) != dims(m2, :lon) || dims(m1, :lat) != dims(m2, :lat)
        msg = "To compute area weigehted RMSE, $m1 and $m2 must be defined on the same lon,lat-grid!"
        throw(ArgumentError(msg))
    end
     # set data set to 0 whenever any of both is missing, 
    # i.e. the returned object has one entry for every model now, if m1, m2 have a 3rd model dimension
    masked_m1 = ifelse.(mask .== true, 0, m1)
    masked_m2 = ifelse.(mask .== true, 0, m2)

    squared_diff = @d (masked_m1 .- masked_m2) .^ 2
    areaweights_mat = makeAreaWeightMatrix(parent(dims(m1, :lon)), parent(dims(m1, :lat)); mask)
    weighted_vals = @d areaweights_mat .* squared_diff
    return sqrt.(sum(coalesce.(weighted_vals, 0), dims=(1,2))[lon=1, lat=1])
end


"""
    normalize(data::Dict{String, <:Number})

Normalize values for every entry in `data` such that they sum up to 1. If remove_zero is
true (default), the returned dictionary does not contain entries for which values were 0.
"""
function normalize(data::Dict{String, <:Number}; remove_zero::Bool=true)
    result = Dict{String, Float64}()
    total = sum(values(data))
    data = remove_zero ? filter(((k, v),) -> v != 0, data) : data
    for k in keys(data)
        result[k] = data[k] ./ total 
    end
    return result
end


"""
    normalizeToYAX(data::Dict{String, <:Number})

Normalize values for every entry in `data` such that they sum up to 1 and return a YAXArray 
with dimension 'diagnostic' whose lookup names are the keys of `data`.
"""
function normalizeToYAX(data::Dict{String, <:Number})
    normed = normalize(data)
    normalized_yax = YAXArray(
        (Dim{:diagnostic}(collect(keys(normed))),), Array{Float64}(undef, length(normed))
    )
    for key in keys(normed)
        normalized_yax[diagnostic = At(key)] = normed[key]
    end
    return normalized_yax
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
    data::AbstractArray; w::Union{YAXArray, Nothing} = nothing, quantiles = [0.167, 0.833]
)
    timesteps = dims(data, :time)
    meta = deepcopy(data.properties)
    meta["_quantiles"] = string.(quantiles)
    uncertainty_ranges = isnothing(w) ?
        YAXArray(
            (dims(data, :time), Dim{:confidence}(["lower", "upper"])),
            Array{Union{Missing, Float64}}(undef, length(timesteps), 2),
            meta
        ) :
        YAXArray(
            (dims(data, :time), dims(w, :weight), Dim{:confidence}(["lower", "upper"])),
            Array{Union{Missing, Float64}}(undef, length(timesteps), length(w.weight), 2),
            meta
        )
    for t in timesteps
        arr = Array(data[time = At(t)])
        if sum(ismissing.(arr)) >= length(arr) - 1 # at least two values must be given
            lower, upper = missing, missing
        else
            if isnothing(w)
                lower, upper = interpolatedWeightedQuantiles(
                    quantiles, collect(skipmissing(arr))
                )
                uncertainty_ranges[time=At(t), confidence=At("lower")] = lower
                uncertainty_ranges[time=At(t), confidence=At("upper")] = upper
            else
                for weight_id in collect(w.weight)
                    weights, df = alignWeightsAndData(data[time = At(t)], w[weight = At(weight_id)])
                    lower, upper = interpolatedWeightedQuantiles(
                        quantiles, collect(skipmissing(df)); weights
                    )
                    uncertainty_ranges[time=At(t), confidence=At("lower"), weight=At(weight_id)] = lower
                    uncertainty_ranges[time=At(t), confidence=At("upper"), weight=At(weight_id)] = upper
                end
            end
        end
    end
    return uncertainty_ranges
end


"""
    interpolatedWeightedQuantiles(
        quantiles::Vector{<:Number}, vals::Vector; weights=nothing
    )

This implementation follows the one used by Brunner et al.
"""
function interpolatedWeightedQuantiles(
    quantiles::Vector{<:Number}, vals::Vector; weights = nothing
)
    if isnothing(weights)
        weights = Array{eltype(vals)}(undef, length(vals))
        weights[ismissing.(vals) .== false] .= 1
    end
    indices_sorted = Array(sortperm(vals)) # gives indices of smallest to largest data point
    weights_sorted = weights[indices_sorted]
    weighted_quantiles = cumsum(weights_sorted) - 0.5 * weights_sorted
    weighted_quantiles = reshape(weighted_quantiles, length(weighted_quantiles), 1)
    # TODO: recheck missing values!
    weighted_quantiles =
        (weighted_quantiles .- minimum(skipmissing(weighted_quantiles))) ./
        maximum(skipmissing(weighted_quantiles))

    interp_linear = Interpolations.linear_interpolation(
        vec(weighted_quantiles),
        vals[indices_sorted],
        extrapolation_bc = Interpolations.Line(),
    )
    return interp_linear(quantiles)
end


"""
    metaDataChecksCMIP(meta::Dict, path::String)

Check model names as retrieved from the metadata for potential inconsistencies wrt filename.
"""
function metaDataChecksCMIP(meta::NCDatasets.CommonDataModel.Attributes, path::String)
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
    dim_symbol = Data.modelDim(data)
    models_data = collect(dims(data, dim_symbol))
    models_weights = collect(dims(weights, dim_symbol))
    data_no_weights = [model for model in models_data if !(model in models_weights)]
    if !isempty(data_no_weights)
        msg = "No weights were computed for follwoing models, thus not considered in the weighted average:"
        @warn msg data_no_weights
        # Only include data for which there are weights
        indices = findall(x -> !(x in data_no_weights), models_data)
        data = Data.indexModel(data, (dim_symbol,), indices)
    end
    # weights for which we don't have data
    weights_no_data = filter(x -> !(x in models_data), models_weights)
    if !isempty(weights_no_data)
        @warn "Weights were renormalized since data of models missing for which weights have been computed: $weights_no_data"
        # renormalize weights
        indices = findall(m -> m in dims(data, dim_symbol), models_weights)
        weights = Data.indexModel(weights, (dim_symbol,), indices)
        weights = weights ./ sum(weights)
    end
    return (weights, data)
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