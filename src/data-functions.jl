"""
    subsetModelData(data::YAXArray, shared_models::Vector{String})

Return subset of `data` containing only data from models in `shared_models`. 

Takes care of metadata.

# Arguments:
- `data`: must have a dimension that contains 'member' or 'model'
- `shared_models`: models, which can either be on level of models or members of models 
('modelname#memberID[_grid]').
"""
function subsetModelData(data::YAXArray, shared_models::Vector{String})
    if isempty(shared_models)
        @warn "Vector of models to subset data to is empty!"
        return data
    end
    data = deepcopy(data)
    level = "model"
    model_dims = modelDims(data)
    dimensions = Array(dims(data, model_dims[1]))
    if level == "member"
        models = map(x -> String(split(x, MODEL_MEMBER_DELIM)[1]), shared_models)
        # if shared_models is on the level of models, the following should be empty
        # otherwise, nothing is filtered out, and members is the same as shared_models 
        members = filter(x -> !(x in models), shared_models)
        if !isempty(members) # shared models on level of members
            indices = findall(m -> m in members, dimensions)
        else
            # important not to use dimensions here, since e.g. model=AWI would be found in dim_names where model is actually AWI-X for instance
            models_data = modelsFromMemberIDs(dimensions)
            indices = findall(m -> m in models, models_data)
        end
    else
        indices = findall(m -> m in shared_models, dimensions)
    end
    data = indexModel(data, model_dims, indices)
    return data
end

"""
    subsetModelData(datamap::DataMap, level::Level=MEMBER_LEVEL)

For those datasets in `datamap` that specify data on the level `level` (i.e. have dimension 
:member or :model), return a new DataMap with subset of data s.t. the new datasets all have 
the same models (MODEL_LEVEL) or members (MEMBER_LEVEL).

If no models are shared across datasets, return the input `datamap`.
"""
function subsetModelData(datamap::DataMap, level::Level = MEMBER_LEVEL)
    shared_models = sharedModels(datamap, level)
    if isempty(shared_models)
        @warn "no shared models in datamap"
        return datamap
    end
    subset = DataMap()
    for (id, data) in datamap
        subset[id] = subsetModelData(deepcopy(data), shared_models)
    end
    return subset
end

function subsetModelData(datamap::DataMap, level::Union{String, Symbol})
    subsetModelData(datamap, getLevel(level))
end



function sharedModels(data::DataMap, level::Level)
    return level == MEMBER_LEVEL ? sharedLevelMembers(data) : sharedLevelModels(data)
end

function sharedModels(data::DataMap, level::Union{String, Symbol})
    return sharedModels(data, getLevel(level))
end

"""
    sharedModels(
        all_paths::Vector{Vector{String}},
        level_shared::Level,
        fn_format::Union{Symbol, String}
    )

# Arguments:
- `all_paths`: every entry refers to the paths to data files for the respective dataset
"""
function sharedModels(
    all_paths::Vector{Vector{String}}, 
    level_shared::Level, 
    fn_format::Union{Symbol, String}
)
    filenames_all_ds = map(paths -> first.(splitext.(basename.(paths))), all_paths)
    filenames_meta_all_ds = map(names -> parseFilename.(names, fn_format), filenames_all_ds)
    models_all = map(meta_data -> map(x -> x.model, meta_data), filenames_meta_all_ds)
    if level_shared == MEMBER_LEVEL
        variants_all = map(meta_data -> map(x -> x.variant, meta_data), filenames_meta_all_ds)
        grids_all = map(meta_data -> map(x -> x.grid, meta_data), filenames_meta_all_ds)
        models_all = map(models_all, variants_all, grids_all) do models, variants, grids
            map(models, variants, grids) do m, v, g
                member = join([m, v], MODEL_MEMBER_DELIM)
                member = !isnothing(g) ? join([member, g], "_") : member
            end
        end
    end
    return reduce(intersect, models_all)
end


"""
    alignPhysics(data::YAXArray, members::Vector{String})

Return new YAXArray with the models retained from `data` that share the same physics as 
the respective model's members in `members`.
All other models that are in `data` but for which no member is specified in `members` are 
also retained.
"""
function alignPhysics(data::YAXArray, members::Vector{String})
    members_data =  Array(dims(data, :member))
    members_kept = Vector{String}()
    for model in modelsFromMemberIDs(members_data; uniq = true)
        allowed_members = filter(m -> startswith(m, model * MODEL_MEMBER_DELIM), members)
        allowed_physics = physicsFromMember.(allowed_members)
        
        members_data_model = filter(x -> startswith(x, model * MODEL_MEMBER_DELIM), members_data)
        physics = physicsFromMember.(members_data_model)
        indices_keep = map(x -> x in allowed_physics, physics)
        if sum(indices_keep) != 0
            push!(members_kept, members_data_model[indices_keep]...)
        end
    end
    return subsetModelData(data, members_kept)

    # if !isnothing(level_shared)
    #     shared_models = sharedModels(data, level_shared, fn_format)
    #     for (id, model_data) in data
    #         data[id] = subsetModelData(model_data, shared_models)
    #     end
    # end
    #return data
end


""" 
    summarizeMembers(data::YAXArray; fn::Function=Statistics.mean)

For each model compute a summary statistic (default: mean) across all its members. 

The returned YAXArray has dimension 'model'.

# Arguments:
- `data::YAXArray`: YAXArray; must have dimension 'member'
- `updateMeta::Bool`: set true if the vectors in the metadata refer to 
different models. Set to false if vectors refer to different variables.
- `fn::Function`: Function to be applied on data
"""
function summarizeMembers(data::YAXArray; fn::Function = Statistics.mean)
    return hasdim(data, :member) ? summarizeMembersVector(data; fn) :
        (hasdim(data, :member1) && hasdim(data, :member2) ? 
        summarizeMembersMatrix(data, true; fn) :
        throw(ArgumentError("Data must have dimension :member or :member1 and :member2. Found: $(dims(data))"))
        )
end

"""
    summarizeMembers!(data::DataMap; fn::Function=Statistics.mean)

Set values for every dataset in `data` to the average across all members of 
each model.
"""
function summarizeMembers!(data::DataMap; fn::Function=Statistics.mean)
    for (k, ds) in data
        if hasdim(ds, :member)
            @info "summarize ensemble members for $k"
            data[k] = summarizeMembers(ds; fn)
        end
    end
    return nothing
end

"""
    summarizeMembers(data::DataMap; fn::Function=Statistics.mean)

Return new DataMap containing every dataset in `data` summarized by applying `fn` 
(default: mean) to all members of each model.
"""
function summarizeMembers(data::DataMap; fn::Function=Statistics.mean)
    df = DataMap()
    for (k, ds) in data
        if hasdim(ds, :member)
            @info "summarize ensemble members for $k"
            df[k] = summarizeMembers(ds; fn)
        end
    end
    return df
end


""" 
    summarizeMembersVector(data::YAXArray; fn::Function=Statistics.mean)

For each model compute a summary statistic (default: mean) across all its members. 
The returned YAXArray has dimension :model (instead of :member).

# Arguments:
- `data::YAXArray`: YAXArray; must have dimension 'member' and at least one other arbitrary dimension 
- `fn::Function`: Function to be applied on data
"""
function summarizeMembersVector(data::YAXArray; fn::Function = Statistics.mean)
    throwErrorIfDimMissing(data, :member)
    if length(dimNames(data)) < 2
        throw(ArgumentError("Vector to summarize members must have at least one other dimension than :member!"))
    end
    data = setLookupsFromMemberToModel(data, ["member"])
    models = Array(dims(data, :model))
    models_uniq = unique(models)
    n_models = length(models_uniq)
    # save indices for each model
    model_indices = Dict{String, Vector{Int}}()
    for (i, m) in pairs(models)
        get!(model_indices, m, Vector{Int}()) 
        push!(model_indices[m], i)
    end
    summarized_data_all = Vector{YAXArray}(undef, n_models)
    for (i, m) in enumerate(models_uniq)
        dat = data[model = model_indices[m]]
        summarized = fn(dat; dims = (:model,))[model = At("combined")]
        meta = subsetMeta(data.properties, model_indices[m]; simplify = true)
        summarized_data_all[i] = YAXArray(dims(summarized), summarized.data, meta)
    end
    summarized_data = combineModelsFromMultipleFiles(
        summarized_data_all; model_names = models_uniq
    )
    summarized_data = replace(summarized_data, NaN => missing)
    return summarized_data
end


"""
    summarizeMembersMatrix(data::AbstractArray, updateMeta::Bool; fn::Function=Statistics.mean)

Compute the average across all members of each model for each given variable 
for model to model data, e.g. distances between model pairs.

# Arguments:
- `data`: with at least dimensions 'member1', 'member2'
- `updateMeta`: set true if the vectors in the metadata refer to different models. 
Set to false if vectors refer to different variables for instance. 
"""
function summarizeMembersMatrix(data::YAXArray, updateMeta::Bool; fn::Function=Statistics.mean)
    throwErrorIfDimMissing(data, [:member1, :member2])
    data = setLookupsFromMemberToModel(data, ["member1", "member2"])
    models = String.(collect(unique(dims(data, :model1))))

    grouped = groupby(data, :model2 => identity)
    averages = map(entry -> mapslices(fn, entry, dims = (:model2,)), grouped) #slow
    combined = cat(averages..., dims = (Dim{:model2}(models)))

    grouped = groupby(combined, :model1 => identity)
    averages = map(entry -> mapslices(fn, entry, dims = (:model1,)), grouped) #slow
    combined = cat(averages..., dims = (Dim{:model1}(models)))

    for m in models
        combined[model1=At(m), model2=At(m)] .= 0
    end

    meta = updateMeta ? updateGroupedDataMetadata(data.properties, grouped) : data.properties
    combined = rebuild(combined; metadata = meta)

    l = Lookups.Categorical(sort(models); order = Lookups.ForwardOrdered())
    combined = combined[model1=At(sort(models)), model2=At(sort(models))]
    combined = DimensionalData.Lookups.set(combined, model1 = l, model2 = l)
    return combined
end



function addMasks!(datamap::DataMap, id_orog_data::String)
    orog_data = datamap[id_orog_data]
    datamap["mask_land"] = getMask(orog_data; mask_out_land = false)
    datamap["mask_ocean"] = getMask(orog_data; mask_out_land = true)
    return nothing
end

### ----------------------------------------------------------------------------------------
###                                LOADING DATA                                           
### ----------------------------------------------------------------------------------------
"""
    loadPreprocData(
        paths::Vector{String},
        filename_format::Union{Symbol, String};
        sorted::Bool = true, 
        dtype::String = "undef",
        names::Vector{String} = Vector{String}(),
        meta_info::Union{Dict{String, T}, Nothing} = nothing,
        constraint_ts::Union{Dict, Nothing} = nothing
    ) where T <: Any

Return data loaded from `paths` as single YAXArray. 

Each path points to a different model, i.e. the data for one dataset is loaded from multiple 
files and all datasets in `paths` must share the same dimensions. Which variable is loaded 
is inferred from the filenames (from `paths`), or if `meta_info` has key 'variable', the 
respective value is used.
"""
function loadPreprocData(
    paths::Vector{String},
    filename_format::Union{Symbol, String};
    sorted::Bool = true, 
    dtype::String = "undef",
    model_names::Vector{String} = Vector{String}(),
    meta_info::Union{Dict{String, T}, Nothing} = nothing,
    constraint_ts::Union{Dict, Nothing} = nothing
) where T <: Any
    is_cmip = lowercase(dtype) == "cmip"
    data = Vector{AbstractArray}(undef, length(paths))
    if isempty(model_names) && is_cmip
        model_names = Vector{String}(undef, length(paths))
    end
    new_dim = is_cmip ? :member : :model
    filenames = first.(splitext.(basename.(paths)))
    filenames_meta = parseFilename.(filenames, filename_format)
    
    for (i, (path, meta)) in enumerate(zip(paths, filenames_meta))
        @debug "processing file.." * path
        ds = NCDataset(path)
        clim_var = !isnothing(meta_info) ? get(meta_info, "variable", meta.variable) : meta.variable
        ds_var = nothing
        try 
            ds_var = ds[clim_var]
        catch e
            @error "Variable '$clim_var' according to filename format '$(string(filename_format))' not found in $path !"
            throw(e)
        end
        # if clim_var == "amoc"
        #     ds_var = ds["msftmz"]
        # end
        props = Dict{String, Any}(ds_var.attrib) # metadata just for this file!
        if clim_var == "msftmz"
            sector = get(ds, "sector", nothing)
            if !isnothing(sector)
                merge!(props, sector.attrib)
            end
        end
        props["path"] = path 
        if is_cmip
            metaDataChecksCMIP(ds.attrib, path)
            mip = get(ds.attrib, "project_id", nothing) # for CMIP5
            mip = isnothing(mip) ? get(ds.attrib, "mip_era", nothing) : mip # for CMIP6
            if isnothing(mip)
                msg = "dtype set to $(dtype). Only CMIP5+CMIP6 supported, yet neither project_id nor mip_era found in data at $path !"
                @warn msg
            else
                member_id_temp = memberIDFromFilenameMeta(meta, mip)
                model_id_temp = String(split(member_id_temp, MODEL_MEMBER_DELIM)[1])
                model_id = fixModelNameMetadata(model_id_temp)
                model_names[i] = replace(member_id_temp, model_id_temp => model_id)
                props["mip_era"] = mip
            end
        end
        dimension_names = dimnames(ds_var)
        dimensions = Vector(undef, length(dimension_names))
        for (idx_dim, d) in enumerate(dimension_names)
            if d in ["bnds", "string21", "time"]
                continue
            end
            dimensions[idx_dim] = Dim{Symbol(d)}(collect(ds_var[d][:]))
        end
        
        exclude_file = false
        if "time" in dimension_names
            # NOTE: just YEAR and MONTH are saved in the time dimension
            times = map(x -> DateTime(Dates.year(x), Dates.month(x)), ds_var["time"][:])
            if absent(constraint_ts) 
                @debug "There were no start and end year for timeseries provided!"
                constraint_ts = Dict()
            end
            indices_time = indicesTimeseries(times, constraint_ts)
            if isempty(indices_time)
                exclude_file = true
            else
                idx_time = findfirst(dimension_names .== "time")
                indices = Vector(undef, length(dimension_names))
                for i in 1:length(dimension_names)
                    indices[i] = i == idx_time ? indices_time : Colon()
                end
                ds_var = ds_var[indices...]
                dimensions[idx_time] = Dim{:time}(collect(times[indices_time]))
            end
        end
        #data[i] = YAXArray(Tuple(dimensions), Array(ds_var), props)
        props["handles"] = ds
        data[i] = exclude_file ? [] : YAXArray(Tuple(dimensions), ds_var, props)
        # replace missing values by NaN?
        # data[i] = YAXArray(Tuple(dimensions), coalesce.(Array(ds_var), NaN), props)
        #close(ds)
    end
    indices = findall(x -> !isempty(x), data)
    if !isempty(model_names)
        model_names = model_names[indices]
    end
    return isempty(indices) ? 
        nothing : 
        combineModelsFromMultipleFiles(data[indices]; model_names, new_dim, sorted, meta=meta_info)
end


"""
    loadDataMapCore(
        all_paths::Vector{Vector{String}},
        ids::Vector{String},
        constraints::Vector{Dict};
        dtype::String = "undef",
        filename_format::Union{Symbol, String} = :cmip,
        sorted::Bool = true,
        preview::Bool = false,
        meta_data::Union{Vector{Dict{String, T}}, Nothing} = nothing
    ) where T <: Any

Load a DataMap instance with keys `ids` that map to data at `all_paths`, where every 
subvector refers the the paths of the data of a single dataset. 

"""
function loadDataMapCore(
    all_paths::Vector{Vector{String}},
    ids::Vector{String},
    constraints::Vector{<:Dict{<:Any, <:Any}};
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :cmip,
    sorted::Bool = true,
    preview::Bool = false,
    meta_data::Union{Vector{<:Dict{String, <:Any}}, Nothing} = nothing
)
    checkConstraint.(constraints)
    if !absent(meta_data) && (length(all_paths) != length(meta_data))
        throw(ArgumentError("size of paths vector and meta data must be equal. Found: paths: $(length(all_paths)), meta_data: $(length(meta_data))"))
    end
    if length(all_paths) != length(ids)
        throw(ArgumentError("size of paths vector and ids must be equal. Found: paths: $(length(all_paths)), ids: $(length(ids))"))
    end
    all_meta = absent(meta_data) ? fill(nothing, length(all_paths)) : meta_data
    data = map(all_paths, all_meta, constraints)  do paths, meta, constraint
        # if provided, do the filtering
        if isnothing(constraint)
            mask = fill(true, length(paths))
        else
            mask = maskFileConstraints(paths, filename_format, constraint)
        end
        constraint_ts = isnothing(constraint) ? nothing : get(constraint, "timeseries", nothing)
        preview ? paths[mask] :
            loadPreprocData(
                paths[mask], filename_format; sorted, dtype, constraint_ts, meta_info = meta 
            )
    end
    if preview
        loaded_data = data
    else
        found_data = map(!isnothing, data)
        if any(x -> x==0, found_data)
            @warn "no data found for ids: $(ids[.!found_data])"
            filter!(!isnothing, data)
            data = YAXArray.(data)
        end
        loaded_data = isempty(data) ? nothing : defineDataMap(data, ids[found_data])
    end
    return loaded_data
end

function loadDataMapCore(
    all_paths::Vector{Vector{String}},
    ids::Vector{String},
    constraint::Union{<:Dict{String, <:Any}, Nothing};
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :cmip,
    sorted::Bool = true,
    preview::Bool = false,
    meta_data::Union{Vector{<:Dict{String, <:Any}}, Nothing} = nothing
)
    constraint = isnothing(constraint) ? Dict{String, Any}() : constraint
    constraints = repeat([constraint], length(all_paths))
    return loadDataMapCore(
        all_paths, ids, constraints; dtype, filename_format, sorted, preview, meta_data
    )
end



"""
    loadDataFromESMValToolRecipes(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool = true,
        constraint::Union{Dict, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true
    )

Load the data from the ESMValTool recipes at `path_recipes` or, if `preview` is true, load
meta data only.

# Arguments:
- `path_data`: top level directory were data is stored.
- `path_recipes`: directory were ESMValTool recipes are stored; these must be the versions 
in the run folder generated by ESMValTool named  RECIPENAME_filled.yml.
- `dir_per_var`: set to true (default) if there is a subdirectory in `path_data` for every 
climate variable to be loaded.
- `constraint`: TODO!
- `preview`: set to true if only meta data should be loaded (default: false).
- `sorted`: if true (default), the data is sorted alphabetically wrt model names.
"""
function loadDataFromESMValToolRecipes(
    path_data::String,
    path_recipes::String;
    dir_per_var::Bool = true,
    constraint::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true, 
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :esmvaltool
)
    checkDataStructure(path_data, dir_per_var)
    meta_data = metaDataFromESMValToolRecipes(path_recipes; constraint)
    paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint)
    if !isnothing(constraint) && !isnothing(get(constraint, "level_shared", nothing))        
        paths = filterPathsSharedModels(paths, constraint["level_shared"], filename_format)
    end
    return loadDataMapCore(
        paths,
        map(x -> x.id, meta_data),
        constraint;
        dtype,
        filename_format,
        sorted,
        preview, 
        meta_data = metadataToDict.(meta_data)
    )
end


"""
    loadDataFromYAML(
        content::Dict;
        constraint::Union{Dict, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true,
        dtype::String = "undef"
    )

Return a DataMap-instance that contains the data specified in `content`, potentially 
constraint by values in `constraint`.

# Arguments:
- `preview::Bool`: if true (default: false), return metadata and corresponding paths without 
actually loading any data.
- `sorted::Bool`: if true (default), model dimension is sorted alphabetically.
- `dtype::String`: if set to "cmip", model dimension of returned data have model names as values.
"""
function loadDataFromYAML(
    yaml_content::Dict;
    constraint::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true, 
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :esmvaltool
)
    fn_err(x) = throw(ArgumentError("$(x) must be provided in config yaml file!"))
    datasets = get(() -> fn_err("datasets"), yaml_content, "datasets")
    base_path = get(() -> fn_err("base_path_data"), yaml_content, "base_path_data")

    all_meta = Vector{}(undef, length(datasets))
    all_paths = Vector{}(undef, length(datasets))
    all_constraints = Vector{}(undef, length(datasets))
    for (i, ds) in enumerate(datasets)
        dir_per_var = get(ds, "dir_per_var", true)
        path_data = joinpath(base_path, get(() -> fn_err("base_dir"), ds, "base_dir"))
        checkDataStructure(path_data, dir_per_var)
        # merge; constraint has precedence constraints defined for individual datasets in the config file
        ds_constraint = get(ds, "subset", Dict())
        if !isnothing(constraint)
            warn_msg = "identical subset keys: arg constraint has precedence over constraint in yaml file!"
            ds_constraint = joinDicts(ds_constraint, constraint; warn_msg) 
        end
        meta_data = metaDataFromYAML(ds)
        paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint=ds_constraint)
        if !isnothing(get(ds_constraint, "level_shared", nothing))        
            msg = "'level_shared' in constraint argument must be one of: $(keys(LEVEL_LOOKUP)), found: $(ds_constraint["level_shared"])."
            level = get(() -> throw(ArgumentError(msg)), LEVEL_LOOKUP, ds_constraint["level_shared"])
            paths = filterPathsSharedModels(paths, level, filename_format)
        end
        all_paths[i] = paths 
        all_meta[i] = meta_data
        # paths and meta_data are vectors! In one loop, several datasets can be loaded (for different variables)
        all_constraints[i] = repeat([ds_constraint], length(paths))
    end
    paths = vcat(all_paths...)
    meta_data = vcat(all_meta...)
    constraints = vcat(all_constraints...)

    # if argument provided, apply subset shared across all datasets
    level_all_data = isnothing(constraint) ? nothing : get(constraint, "level_shared", nothing)
    if !isnothing(level_all_data)
        paths = filterPathsSharedModels(paths, level_all_data, filename_format)
    end
    return loadDataMapCore(
        paths, 
        map(x -> x.id, meta_data),
        constraints; 
        dtype, 
        filename_format, 
        sorted, 
        preview, 
        meta_data = metadataToDict.(meta_data)
    )
end

function defineDataMap(
    yaml_content::Dict;
    constraint::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true,
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :esmvaltool
)
    return loadDataFromYAML(yaml_content; constraint, preview, sorted, dtype, filename_format)
end


function defineDataMap(
    path_config::String;
    constraint::Union{Dict,Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true,
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :esmvaltool
)
    return loadDataFromYAML(
        YAML.load_file(path_config); constraint, preview, sorted, dtype, filename_format
    )
end


function defineDataMap(
    path_data::String,
    path_recipes::String,
    source::Symbol;
    dir_per_var::Bool = true,
    constraint::Union{Dict, Nothing} = nothing,
    preview::Bool = false,
    sorted::Bool = true, 
    dtype::String = "undef",
    filename_format::Union{Symbol, String} = :esmvaltool
)
    if source != :esmvaltool_recipes
        throw(ArgumentError("To load data from esmvaltool recipes, set source argument to :esmvaltool_recipes, found: $(source)."))
    end
    return loadDataFromESMValToolRecipes(
        path_data, path_recipes; dir_per_var, constraint, preview, sorted, dtype, filename_format
    )
end


# Loading data directly from given directories
"""
    defineDataMap(
        path_data_dir::String, 
        id::String;
        meta_data::Dict{String, T} = Dict{String, Any}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::String = "undef",
        filename_format::Union{Symbol, String} = :cmip
    ) where T <: Any

Return DataMap with entry `id` with the data (model and/or obs depending on `dtype`) from 
all .nc files in `paths` and all .nc files in all directories in `paths`, possibly 
constraint by `constraint`.
"""
function defineDataMap(
    paths::Vector{String}, 
    id::String;
    meta_data::Dict{String, T} = Dict{String, Any}(),
    constraint::Union{Dict, Nothing} = nothing,
    filename_format::Union{Symbol, String} = :cmip,
    dtype::String = "undef",
    preview::Bool = false,
    sorted::Bool = true
) where T <: Any
    paths_dirs = filter(x -> isdir(x), paths)
    paths_ncfiles = filter(x -> isfile(x) && endswith(x, ".nc"), paths)
    paths_to_files = vcat(collectNCFilePaths.(paths_dirs)..., paths_ncfiles)
    return loadDataMapCore(
        [paths_to_files], 
        [id],
        constraint;
        dtype, 
        filename_format,
        sorted,
        preview,
        meta_data = [meta_data]
    ) 
end


"""
    defineDataMap(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::String = "undef",
        filename_format::Union{Symbol, String} = :cmip
    ) where T <: Any

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_data_dirs`, possibly constraint by `constraint`.
"""
function defineDataMap(
    paths_data_dirs::Vector{String}, 
    data_ids::Vector{String};
    meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
    constraint::Union{Dict, Nothing} = nothing,
    filename_format::Union{Symbol, String} = :cmip,
    dtype::String = "undef",
    preview::Bool = false,
    sorted::Bool = true,
) where T <: Any
    paths_to_files = collectNCFilePaths.(paths_data_dirs)
    if !isnothing(constraint) && !isempty(constraint) && !isnothing(get(constraint, "level_shared", nothing))        
        paths_to_files = filterPathsSharedModels(
            paths_to_files, constraint["level_shared"], filename_format
        )
    end
    return loadDataMapCore(
        paths_to_files, 
        data_ids,
        constraint; 
        dtype, 
        filename_format, 
        sorted, 
        preview,
        meta_data
    )
end


"""
    defineDataMap(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
        constraint::Union{Dict, Nothing} = nothing,
        sorted::Bool = true, 
        dtype::String = "undef",
        filename_format::Union{Symbol, String} = :cmip
    ) where T <: Any

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_dirs`, possibly constraint by `constraint`.

For every loaded dataset (entry in built DataMap), files are loaded from several directories;
each entry in `paths_data_dirs` points the the vector of data directories from where data 
is loaded for that dataset.
"""
function defineDataMap(
    paths_data_dirs::Vector{Vector{String}}, 
    data_ids::Vector{String};
    meta_data::Vector{Dict{String, T}} = Vector{Dict{String, Any}}(),
    constraint::Union{Dict, Nothing} = nothing,
    filename_format::Union{Symbol, String} = :cmip,
    dtype::String = "undef",
    preview::Bool = false,
    sorted::Bool = true
) where T <: Any
    paths_to_files = map(paths_data_dirs) do paths 
        vcat(collectNCFilePaths.(paths)...)
    end
    if !isnothing(constraint) && !isempty(constraint) && !isnothing(get(constraint, "level_shared", nothing))        
        paths_to_files = filterPathsSharedModels(
            paths_to_files, constraint["level_shared"], filename_format
        )
    end
    return loadDataMapCore(
        paths_to_files, 
        data_ids,
        constraint; 
        dtype, 
        filename_format, 
        sorted, 
        preview,
        meta_data
    )
end
