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
        @warn "Vector of models to subset data to is empty! Data is not filtered."
        return data
    end
    model_dims = modelDims(data)
    dimensions = Array(dims(data, model_dims[1]))
    if occursin("member", string(model_dims[1]))
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
    subsetModelData(dm::DataMap, level::Symbol=:member; ids::Vector{String}=Vector{String}())

For those datasets in `dm` that specify data on the level `level` (i.e. have dimension 
:member or :model), return a new DataMap with subset of data s.t. the new datasets all have 
the same models or members.

If no models are shared across datasets, return the input `dm`.
"""
function subsetModelData(
    dm::DataMap, level::Symbol = :member; ids::Vector{String}=Vector{String}()
)
    shared_models = sharedModels(dm, toLevel(Val(level)))
    if isempty(shared_models)
        @warn "no shared models in input datamap!"
        return dm
    end
    subset = DataMap()
    if isempty(ids)
        ids = collect(keys(dm))
    end
    for id in ids
        subset[id] = subsetModelData(deepcopy(dm[id]), shared_models)
    end
    return subset
end

# function subsetModelData(dm::DataMap, level::Symbol; ids::Vector{String}=Vector{String}())
#     subsetModelData(dm, getLevel(level); ids)
# end


"""
    Subset all entries in dm to the shared set of models on 'level' among all entries.
"""
function subsetModelData!(
    dm::DataMap, level::Symbol = :member; ids::Vector{String} = Vector{String}()
)
    shared_models = sharedModels(dm, toLevel(Val(level)))
    if isempty(shared_models)
        @warn "no shared models in input datamap!"
        return dm
    end
    if isempty(ids)
        ids = collect(keys(dm))
    end
    for id in ids
        dm[id] = subsetModelData(dm[id], shared_models)
    end
    return nothing
end

# function subsetModelData!(dm::DataMap, level::Symbol; ids::Vector{String}=Vector{String}())
#     subsetModelData!(dm, getLevel(level); ids)
# end


function sharedModels(data::DataMap, level::Level)
    return isa(level, LevMember) ? sharedLevelMembers(data) : sharedLevelModels(data)
end

# function sharedModels(data::DataMap, level::Symbol)
#     return sharedModels(data, getLevel(level))
# end

# """
#     sharedModels(
#         all_paths::Vector{Vector{String}},
#         level_shared::Level,
#         fn_format::Symbol
#     )

# # Arguments:
# - `all_paths`: every entry refers to the paths to data files for the respective dataset
# """
# function sharedModelsFromPaths(
#     all_paths::Vector{Vector{String}}, level::Level, fn_format::FilenameFormat
# )
#     filenames_all_ds = map(paths -> first.(splitext.(basename.(paths))), all_paths)
#     filenames_meta_all_ds = map(names -> parseFilename.(names, fn_format), filenames_all_ds)
#     models_all = map(meta_data -> map(x -> x.model, meta_data), filenames_meta_all_ds)
#     if isa(level, LevMember)
#         variants_all = map(meta_data -> map(x -> x.variant, meta_data), filenames_meta_all_ds)
#         grids_all = map(meta_data -> map(x -> x.grid, meta_data), filenames_meta_all_ds)
#         models_all = map(models_all, variants_all, grids_all) do models, variants, grids
#             map(models, variants, grids) do m, v, g
#                 member = join([m, v], MODEL_MEMBER_DELIM)
#                 member = !isnothing(g) ? join([member, g], "_") : member
#             end
#         end
#     end
#     return reduce(intersect, models_all)
# end

function sharedModelsFromPaths(
    all_paths::Vector{Vector{String}}, level::Level, fn_format::AbstractFnFormat
)
    models_all = Vector{Set{String}}(undef, length(all_paths))
    @inbounds for (i,paths) in pairs(all_paths)
        models = Set{String}()
        @inbounds for p in paths
            #name = first(splitext(basename(p)))
            meta = parseFilename(p, fn_format)
            push!(models, _model_key(meta, level))
        end
        models_all[i] = models
    end
    return collect(reduce(intersect, models_all))
end

# dispatch instead of if/isa branching
_model_key(meta, ::LevModel) = meta.model

function _model_key(meta, ::LevMember)
    key = string(meta.model, MODEL_MEMBER_DELIM, meta.variant)
    return isnothing(meta.grid) ? key : string(key, "_", meta.grid)
end


function sharedModelsFromMeta(meta_data::Vector{Vector{FilenameMeta}}, level::Level)
    models_per_dataset = Vector{Set{String}}(undef, length(meta_data))
    @inbounds for (i, meta_data_files) in pairs(meta_data)
        models = Set{String}()
        @inbounds for meta in meta_data_files            
            push!(models, _model_key(meta, level))
        end
        models_per_dataset[i] = models
    end
    return collect(reduce(intersect, models_per_dataset))
end




"""
    alignPhysics(data::YAXArray, members::Vector{String})

Return new YAXArray with the models retained from `data` that share the same physics as 
the respective model's members in `members`.
All other models that are in `data` but for which no member is specified in `members` are 
also retained.
"""
function alignPhysics(data::YAXArray, members::AbstractVector{String})
    members_data = lookup(data, :member)
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
    summarizeMembers(
        data::AbstractArray,
        dim::Int,
        members::AbstractArray{String}; 
        fn::Function=Statistics.mean
    )

For each model compute a summary statistic (default: mean) across all its members. 
Return Array summarized in model dimension.

# Arguments:
- `data::AbstractArray`: must have dimension 'member' and at least one other arbitrary dimension 
- `dim::Int`: index of member dimension
- `members::AbstractArray{String}`: names of members corresponding to member dimension in `data`.
- `fn::Function`: Function to be applied on data
"""
function summarizeMembers(
    data::AbstractArray, 
    dim::Int, 
    members::AbstractArray{String}; 
    fn::Function = Statistics.mean
)
    models_all = membersToModels(members)
    models_uniq = unique(models_all)
    n_models = length(models_uniq)
    # save indices for each model
    model_indices = Dict{String, Vector{Int}}()
    for (i, m) in pairs(models_all) # pairs maps from indices to entries
        get!(model_indices, m, Vector{Int}()) 
        push!(model_indices[m], i)
    end
    T = eltype(data)
    summarized_data_all = Vector{AbstractArray{T}}(undef, n_models)
    for (i, m) in enumerate(models_uniq)
        dat = selectdim(data, dim, model_indices[m])
        summarized = fn(dat; dims = (dim,))
        summarized_data_all[i] = selectdim(summarized, dim, 1)
    end
    n_dims = length(size(data))
    summarized_data = cat(summarized_data_all..., dims = n_dims)
    if n_dims > 1
        indices = Array(1:n_dims)
        indices[end] = dim
        indices[dim] = n_dims
        summarized_data = permutedims(summarized_data, Tuple(indices))
    end
    return summarized_data 
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
    data = setLookupsFromMemberToModel(data, [:member])
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

    idx_model = indexDim(data, :model)
    for (i, m) in enumerate(models_uniq)
        #dat = data[model = model_indices[m]]
        dat = selectdim(data, idx_model, model_indices[m]) # this is more than 2x better in terms of allocations

        summarized = fn(dat; dims = (:model,))[model = At("combined")]
        #summarized = selectdim(fn(dat; dims = (:model,)), idx_model, 1)

        meta = subsetMeta(deepcopy(data.properties), model_indices[m]; simplify = true)
        summarized_data_all[i] = YAXArray(dims(summarized), summarized.data, meta)
    end
    summarized_data = combineModelsFromMultipleFiles(
        summarized_data_all; model_names = models_uniq
    )
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
    data = setLookupsFromMemberToModel(data, [:member1, :member2])
    models_all = collect(dims(data, :model1))
    models = unique(models_all)
    other_dims = otherdims(data, (:model1, :model2))
    mat = YAXArray(
        (Dim{:model1}(models), Dim{:model2}(models), other_dims...), 
        zeros(length(models), length(models), length.(other_dims)...),
        deepcopy(data.properties)
    )
    for (i, m1) in enumerate(models[1: end-1])
        for m2 in models[i+1 : end]
            indices1 = findall(x -> x == m1, models_all)
            indices2 = findall(x -> x == m2, models_all)
            slice = data[model1 = indices1, model2 = indices2]
            summarized = fn(slice, dims = (:model1, :model2))[model1=At("combined"), model2=At("combined")]
            
            mat[model1 = At(m1), model2 = At(m2)] = isempty(other_dims) ? Array(summarized)[1] : summarized
            # build up symmetric matrix
            mat[model1 = At(m2), model2 = At(m1)] = isempty(other_dims) ? Array(summarized)[1] : summarized 
        end
    end
    return mat
    # TODO: check metadata update!
    # meta = updateMeta ? updateGroupedDataMetadata(data.properties, grouped) : data.properties
    # combined = rebuild(combined; metadata = meta)
    # l = Lookups.Categorical(sort(models); order = Lookups.ForwardOrdered())
    # combined = combined[model1=At(sort(models)), model2=At(sort(models))]
    # combined = DimensionalData.Lookups.set(combined, model1 = l, model2 = l)
    # return combined
end


function addMasks!(datamap::DataMap, id_orog_data::String)
    orog_data = datamap[id_orog_data]
    datamap["mask_land"] = getMask(orog_data; mask_out_land = false)
    datamap["mask_ocean"] = getMask(orog_data; mask_out_land = true)
    return nothing
end


function filterOutModels(data::YAXArray, excluded_models::Vector{String})
    throwErrorIfDimMissing(data, [:model, :member]; include = :any)
    model_dim = modelDim(data)
    models = lookup(data, model_dim)
    indices_ok = findall(
        x -> !(x in excluded_models) && !(split(x, MODEL_MEMBER_DELIM)[1] in excluded_models), models
    )
    return model_dim == :model ? data[model = indices_ok] : data[member = indices_ok]
end


function listModels(data::YAXArray)
    model_dim = modelDim(data)
    if model_dim == :member
        models = modelsFromMemberIDs(data; uniq=false)
        members = membersFromMemberIDs(data)
    else 
        models = lookup(data, model_dim)
        members = fill("", length(models))
    end
    return DataFrame(model = collect(models), run = members)
end

### ----------------------------------------------------------------------------------------
###                                LOADING DATA                                           
### ----------------------------------------------------------------------------------------
"""
    loadPreprocData(
        paths::Vector{String},
        filename_format::Symbol;
        sorted::Bool = true, 
        dtype::String = "cmip",
        model_names::Vector{String} = Vector{String}(),
        meta_info::Union{Dict{String, String}, Nothing} = nothing,
        constraint_ts::Dict{String, Int} = Dict{String, Int}()
    )

Return data loaded from `paths` as single YAXArray. 

Each path points to a different model, i.e. the data for one dataset is loaded from multiple 
files and all datasets in `paths` must share the same dimensions. Which variable is loaded 
is inferred from the filenames (in `paths`), or if `meta_info` has key 'variable', the 
respective value is used.
"""
function loadPreprocData(
    meta_data::Vector{FilenameMeta};
    #filename_format::AbstractFnFormat;
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}} = (start_year = typemin(Int), end_year = typemax(Int)),
    dtype::String = "cmip",
    sorted::Bool = true,
    meta_info::Union{Dict{String, String}, Nothing} = nothing
)
    data = YAXArray[]
    model_names = String[]
    new_dim = dtype == "cmip" ? :member : :model
    
    # iterate over meta data for each file
    @inbounds for (i, meta) in enumerate(meta_data)
        @debug "processing file $(meta.path) ..."
        ds = open_dataset(meta.path)
        #clim_var = !isnothing(meta_info) ? get(meta_info, "variable", meta.variable) : meta.variable
        ds_var = try
            ds[Symbol(meta.variable)]
        catch e
            if e isa KeyError
                throw(ArgumentError("Variable $(meta.variable) not found in $(meta.path)!"))
            else
                rethrow()
            end
        end
        # if clim_var == "amoc"
        #     ds_var = ds["msftmz"]
        # end        
        # apply timeseries constraint
        dimension_names = dimNames(ds_var)

        ndim = length(dimension_names)
        dimensions = Vector(undef, ndim)
        for (idx_dim, d) in enumerate(dimension_names)
            if d != :time
            # if d in ["bnds", "string21", "time"]
            #     continue
            # end
                dimensions[idx_dim] = Dim{d}(Array(dims(ds_var, d)))
                #dimensions[idx_dim] = Dim{Symbol(d)}(collect(ds_var[d][:]))
            end
        end

        exclude_file = false
        if :time in dimension_names
            # NOTE: just YEAR is saved in the time dimension
            times = [DateTime(Dates.year(x), Dates.month(x)) for x in lookup(ds_var, :time)]
            #add meta data for time (necessary?)
            # time_meta = Dict{String, Any}()
            # props["time-meta"] = time_meta
            # for t in (:year, :month)
            #     if hasproperty(ds, t)
            #         var = ds[t]
            #         for (k,v) in var.properties
            #             props["time-meta"][k] = v
            #         end
            #     end
            # end
            indices_time = indicesTimeseries(times, constraint_ts)
            if isempty(indices_time)
                exclude_file = true
            else
                idx_time = findfirst(==( :time ), dimension_names)
                indices = ntuple(i -> i == idx_time ? indices_time : Colon(), length(dimension_names))
                ds_var = ds_var[indices...]
                dimensions[idx_time] = Dim{:time}(collect(times[indices_time]))
            end
        end
        #TODO: This is apparently not used!!!!
        # fv = get(props, "_FillValue", missing)
        # ds_var = map(x -> x == fv  ? missing : x, ds_var)
        # raw = Array(ds_var)
        # if haskey(props, "_FillValue")
        #     fv = props["_FillValue"]
        #     #raw = Array(ds_var)
        #     Td = eltype(raw)
        #     ds_var_out = Array{Union{Missing, Td}}(undef, size(raw))
        #     @inbounds for i in eachindex(raw)
        #         x = raw[i]
        #         ds_var_out[i] = isequal(x, fv) ? missing : x
        #     end
        #     raw = ds_var_out
        # end                
        if !exclude_file
            props = copy(ds_var.properties) # metadata just for this file!
            if meta.variable == "msftmz"
                sector = get(ds, "sector", nothing)
                if !isnothing(sector)
                    merge!(props, sector.properties)
                end
            end
            props["path"] = meta.path
            props["mip_era"] = meta.mip
            if dtype == "cmip"
                # returns member name
                push!(model_names, fixModelNameInconsistenciesCMIP(ds.properties, meta))
            end
            push!(data, YAXArray(Tuple(dimensions), allowmissing(ds_var), props))
        end
    end
    return isempty(data) ? nothing : 
        combineModelsFromMultipleFiles(data; model_names, new_dim, sorted, meta=meta_info)
end



function checkInput(
    all_paths::Vector{T}, 
    ids::Vector{String};
    meta_data::Vector{Dict{String, String}} = Dict{String, String}[]
) where {T}
    if !absent(meta_data) && (length(all_paths) != length(meta_data))
        throw(ArgumentError("size of paths vector and meta data must be equal. Found: paths: $(length(all_paths)), meta_data: $(length(meta_data))"))
    end
    np, ni = length.([all_paths, ids])
    if np != ni
        throw(ArgumentError("'all_paths' and 'ids' must have the same length, found: ($(length(all_paths)) vs. $(length(ids)))."))
    end
end


function _previewDataMapCore(
    meta_data::Vector{Vector{FilenameMeta}}, ids::AbstractArray{String}
)
    preview_map = PreviewMap()
    for i in eachindex(ids)
        preview_map[ids[i]] = meta_data[i] 
    end
    return preview_map
end

function _loadDataMapCore(
    meta_data_per_dataset::Vector{Vector{FilenameMeta}},
    ids::Vector{String};
    dtype::String = "cmip",
    fn_format::AbstractFnFormat = FF_CMIP(),
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}} = (start_year = typemin(Int), end_year = typemax(Int)),
    level::AbstractLevel = NoLevel(),
    sorted::Bool = true,
    #meta_info::Vector{Dict{String, String}} = Dict{String, String}[]
)
    # TODO: handle meta info
    data = Vector{YAXArray}(undef, length(meta_data_per_dataset))
    indices_found = Int[];
    for (i, meta_data) in enumerate(meta_data_per_dataset)
        df = loadPreprocData(
            meta_data;
            #fn_format;
            constraint_ts, 
            dtype, 
            sorted, 
            #meta_info = meta_info  
        )
        if !isnothing(df)
            data[i] = df
            push!(indices_found, i)
        end
    end
    data = data[indices_found]
    # no_data_found = emptyYAX.(data)
    # data_found = .!no_data_found
    # if any(no_data_found)
    #     #@warn "no data found for ids: $(ids[.!found_data])"
    #     #filter!(!isnothing, data)
    #     data = data[data_found]
    #     #data = YAXArray.(data) # why needed?
    # end
    if isempty(data)
        # no data was found at all
        return nothing # TODO: better to return an empty array?!
    end
    # if isa(level, Level) && dtype == "cmip"
    #     models = isa(level, LevModel) ?  modelsFromMemberIDs.(data; uniq=true) : Array.(lookup.(data, :member))
    #     shared_models = reduce(intersect, models)
    #     data = map(ds -> subsetModelData(ds, shared_models), data)
    # end
    #ids_filtered = ids[data_found]
    ids_filtered = ids[indices_found]
    return defineDataMap(data, ids_filtered)
end

# previous:
# """
#     _loadDataMapCore(
#         all_paths::Vector{Vector{String}},
#         ids::Vector{String},
#         constraint::Union{Dict{String, <:AbstractArray{String}}, Nothing} = nothing,
#         constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#         level::Union{Symbol, Nothing} = nothing,
#         dtype::String = "cmip",
#         filename_format::Symbol = :cmip,
#         sorted::Bool = true,
#         meta_data::Union{Vector{Dict{String, String}}, Nothing} = nothing
#     )

# Load a DataMap instance with keys `ids` that map to data at `all_paths`, where every 
# subvector refers the paths of the data of a single dataset. 


# # Arguments:
# - `all_paths`: each subvector should contain paths to directories that store data for the same experiment and variable.
# """
# function _loadDataMapCore(
#     all_paths::Vector{Vector{String}},
#     ids::Vector{String};
#     constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#     constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#     level::AbstractLevel = LevNone(),
#     dtype::String = "cmip",
#     fn_format::AbstractFnFormat = FF_CMIP(),
#     sorted::Bool = true,
#     meta_data::Vector{Dict{String, String}} = Dict{String, String}[]
# )
#     #checkConstraint.(constraints)
#     #checkInput(all_paths, ids; meta_data)

#     all_meta = absent(meta_data) ? fill(nothing, length(all_paths)) : meta_data
#     #@info "length(all_paths): $(length(all_paths)); length(all_paths[1]): $(length(all_paths[1]))"
#     # TODO: add possibility to have one constraint per dataset 
#     #data = Vector{Vector{YAXArray}}(undef, length(all_paths))
#     data = Vector{YAXArray}(undef, length(all_paths))
#     for (i, (paths, meta)) in enumerate(zip(all_paths, all_meta))
#         data[i] = loadPreprocData(
#             paths, 
#             fn_format; 
#             constraint_ts, 
#             dtype, 
#             sorted, 
#             meta_info = meta 
#         )
#     end
#     no_data_found = emptyYAX.(data)
#     data_found = .!no_data_found
#     if any(no_data_found)
#         #@warn "no data found for ids: $(ids[.!found_data])"
#         #filter!(!isnothing, data)
#         data = data[data_found]
#         #data = YAXArray.(data) # why needed?
#     end
#     if isempty(data)
#         # no data was found at all
#         return nothing # TODO: better to return an empty array?!
#     end
#     # if isa(level, Level) && dtype == "cmip"
#     #     models = isa(level, LevModel) ?  modelsFromMemberIDs.(data; uniq=true) : Array.(lookup.(data, :member))
#     #     shared_models = reduce(intersect, models)
#     #     data = map(ds -> subsetModelData(ds, shared_models), data)
#     # end
#     ids_filtered = ids[data_found]
#     # return defineDataMap(data, ids_filtered)
#     return defineDataMap(data, ids_filtered)
# end

# function loadDataMapCore(
#     all_paths::Vector{Vector{String}},
#     ids::Vector{String},
#     constraint::Union{Dict{String, T}, Nothing};
#     dtype::String = "cmip",
#     filename_format::Symbol = :cmip,
#     sorted::Bool = true,
#     preview::Bool = false,
#     meta_data::Union{Vector{Dict{String, String}}, Nothing} = nothing
# ) where T <: AbstractArray{String}
#     constraint = isnothing(constraint) ? Dict{String, Vector{String}}() : constraint
#     constraints = repeat([constraint], length(all_paths))
#     return loadDataMapCore(
#         all_paths, ids, constraints; dtype, filename_format, sorted, preview, meta_data
#     )
# end

# function _applySharedLevel(
#     all_paths::AbstractVector{<:AbstractVector{String}}, level::Level, fn_format::FilenameFormat
# )
#     return filterPathsSharedModels(all_paths, level, fn_format)
# end

function _applyLevel!(
    meta_data_per_dataset::AbstractVector{<:AbstractVector{FilenameMeta}}, level::Level
)
    shared = sharedModelsFromMeta(meta_data_per_dataset, level)
    indices_members = findall(x -> occursin(MODEL_MEMBER_DELIM, x), shared)
    indices_models = [i for i in 1:length(shared) if !(i in indices_members)]
    #shared_members = shared[indices_members] # TODO
    shared_models = shared[indices_models]
    for (i, meta_files) in enumerate(meta_data_per_dataset)
        # handle members TODO
        #handle models
        mask = [meta.model in shared_models for meta in meta_files]
        meta_data_per_dataset[i] = meta_files[mask]
    end
    return nothing
end

function _applyConstraint!(
    meta_data_per_dataset::AbstractVector{<:AbstractVector{FilenameMeta}},
    constraint::Constraint
)
    for (i, meta_files) in enumerate(meta_data_per_dataset)
        mask = [isRetained(meta, constraint) for meta in meta_files]
        meta_data_per_dataset[i] = meta_files[mask]
    end
    return nothing
end

function _applyConstraintTS!(
    meta_data_per_dataset::AbstractVector{<:AbstractVector{FilenameMeta}},
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}}
)
    start_cs = constraint_ts.start_year
    end_cs = constraint_ts.end_year
    if end_cs < start_cs
        throw(ArgumentError("End of timseries must be after start!"))
    end

    for (i, meta_files) in enumerate(meta_data_per_dataset)
        mask = Vector{Bool}(undef, length(meta_files))
        for (j, meta) in enumerate(meta_files)
            # timerange given as string, e.g. 185001-185012
            start_file = parse(Int, meta.timerange[1:4])
            end_file = parse(Int, meta.timerange[8:11])
            mask[j] = (start_file >= start_cs && start_file <= end_cs) ||
                (end_file >= start_cs && end_file <= end_cs)
        end
        meta_data_per_dataset[i] = meta_files[mask]
    end
    return nothing
end


"""
    loadDataFromESMValToolRecipes(
        path_data::String,
        path_recipes::String;
        dir_per_var::Bool = true,
        constraint::Union{Dict{String, <:AbstractArray{String}}, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true,
        dtype::String = "cmip",
        filename_format::Symbol = :esmvaltool
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
    constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
    constraint_ts::Dict{String, Int} = Dict{String, Int}(),
    level::AbstractLevel = NoLevel(),
    dtype::String = "cmip",
    fn_format::AbstractFnFormat = ESMVTFormat(),
    sorted::Bool = true,
    preview::Bool = false
)
    checkDataStructure(path_data, dir_per_var)
    meta_data = metaDataFromESMValToolRecipes(path_recipes; constraint)
    paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint)
    if !isnothing(level)       
        paths = filterPathsSharedModels(paths, level, fn_format)
    end
    # if preview
    #     return _previewDataMapCore(paths; constraint, level, dtype, filename_format)
    # else
    return _loadDataMapCore(
        paths,
        getfield.(meta_data, :id);
        constraint,
        constraint_ts,
        level,
        dtype,
        fn_format,
        sorted,
        meta_data = metadataToDict.(meta_data)
    )
    # end
end


"""
    loadDataFromYAML(
        yaml_content::Dict;
        constraint::Union{Dict{String, <:AbstractArray{String}}, Nothing} = nothing,
        preview::Bool = false,
        sorted::Bool = true,
        dtype::String = "cmip",
        filename_format::Symbol = :esmvaltool
    )

Return a DataMap-instance that contains the data specified in `yaml_content`, potentially 
constraint by values in `constraint`.

# Arguments:
- `preview::Bool`: if true (default: false), return metadata and corresponding paths without 
actually loading any data.
- `sorted::Bool`: if true (default), model dimension is sorted alphabetically.
- `dtype::String`: if set to "cmip", model dimension of returned data have model names as values.
"""
function loadDataFromYAML(
    yaml_content::Dict;
    constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
    constraint_ts::Dict{String, Int} = Dict{String, Int}(),
    level::AbstractLevel = NoLevel(),
    dtype::String = "cmip",
    fn_format::AbstractFnFormat = ESMVTFormat(),
    sorted::Bool = true, 
    preview::Bool = false
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
        # merge; constraint has precedence over constraints defined for individual datasets in the config file
        ds_constraint = get(ds, "subset", Dict())
        if !isnothing(constraint)
            warn_duplicates = "identical subset keys: arg constraint has precedence over constraint in yaml file!"
            ds_constraint = joinDicts(ds_constraint, constraint; warn_msg = warn_duplicates) 
        end
        # TODO: constraint_ts now just possible via argument, not inside yaml, should be done here also for timeseries constraints

        meta_data = metaDataFromYAML(ds)
        paths = resolvePathsFromMetaData.(meta_data, path_data, dir_per_var; constraint=ds_constraint)
        if isa(level, Level)        
            paths = filterPathsSharedModels(paths, level, fn_format)
        end
        all_paths[i] = paths 
        all_meta[i] = meta_data
        # paths and meta_data are vectors! In one loop, several datasets can be loaded (for different variables)
        #all_constraints[i] = repeat([ds_constraint], length(paths))
        all_constraints[i] = ds_constraint
    end
    meta_data = vcat(all_meta...)

    # if preview 
    #     return _previewDataMapCore(vcat(all_paths...); constraint, level, dtype, filename_format)
    # else 
    return _loadDataMapCore(
        vcat(all_paths...), 
        getfield.(meta_data, :id);
        constraint = constraint, #all_constraints, # TODO: now only defined for single constraint! former: #vcat(all_constraints...), 
        constraint_ts,
        level,
        dtype, 
        fn_format, 
        sorted, 
        meta_data = metadataToDict.(meta_data)
    )
    # end
end

# function defineDataMap(
#     yaml_content::Dict;
#     constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#     constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#     level::Symbol = :none,
#     dtype::String = "cmip",
#     filename_format::Symbol = :esmvaltool,
#     sorted::Bool = true,
#     preview::Bool = false
# )
#     level = toLevel(Val(level))
#     fn_format = toFF(Val(filename_format))
#     return loadDataFromYAML(
#         yaml_content; 
#         constraint, 
#         constraint_ts, 
#         level, 
#         dtype, 
#         fn_format, 
#         sorted, 
#         preview
#     )
# end


# function defineDataMap(
#     path_config::String;
#     constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#     constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#     level::Symbol = :none,
#     dtype::String = "cmip",
#     filename_format::Symbol = :esmvaltool,
#     sorted::Bool = true,
#     preview::Bool = false
# )
#     level = toLevel(Val(level))
#     fn_format = toFF(Val(filename_format))
#     return loadDataFromYAML(
#         YAML.load_file(path_config); 
#         constraint, 
#         constraint_ts,
#         level,
#         dtype, 
#         fn_format, 
#         sorted,
#         preview
#     )
# end


# function defineDataMap(
#     path_data::String,
#     path_recipes::String,
#     source::Symbol;
#     dir_per_var::Bool = true,
#     constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#     constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#     level::Symbol = :none,
#     dtype::String = "cmip",
#     filename_format::Symbol = :esmvaltool,
#     sorted::Bool = true,
#     preview::Bool = false
# )
#     if source != :esmvaltool_recipes
#         throw(ArgumentError("To load data from esmvaltool recipes, set source argument to :esmvaltool_recipes, found: $(source)."))
#     end
#     level = toLevel(Val(level))
#     fn_format = toFF(Val(filename_format))
    
#     return loadDataFromESMValToolRecipes(
#         path_data, 
#         path_recipes; 
#         dir_per_var, 
#         constraint, 
#         constraint_ts, 
#         level, 
#         dtype, 
#         fn_format, 
#         sorted, 
#         preview
#     )
# end


# # Loading data directly from given directories
# """
#     defineDataMap(
#         paths::Vector{String}, 
#         id::String;
#         meta_data::Dict{String, String} = Dict{String, String}(),
#         constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#         filename_format::Symbol = :cmip
#         dtype::String = "cmip",
#         preview::Bool = false,
#         sorted::Bool = true,
#     )

# Return DataMap with entry `id` with the data (model and/or obs depending on `dtype`) from 
# all .nc files in `paths` and all .nc files in all directories in `paths`, possibly 
# constraint by `constraint`.
# """
# function defineDataMap(
#     paths::Vector{String}, 
#     id::String;
#     constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
#     constraint_ts::Dict{String, Int} = Dict{String, Int}(),
#     level::Symbol = :none,
#     dtype::String = "cmip",
#     filename_format::Symbol = :cmip,
#     sorted::Bool = true,
#     meta_data::Dict{String, String} = Dict{String, String}()
# )
#     level = toLevel(Val(level))
#     fn_format = toFF(Val(filename_format))

#     paths_dirs = filter(x -> isdir(x), paths)
#     paths_ncfiles = filter(x -> isfile(x) && endswith(x, ".nc"), paths)
#     paths_to_files = vcat(collectNCFilePaths.(paths_dirs)..., paths_ncfiles)

#     paths_to_files = _applyConstraint(paths_to_files, constraint; level, dtype, fn_format)

#     return _loadDataMapCore(
#         [paths_to_files], 
#         [id];
#         constraint,
#         constraint_ts,
#         level,
#         dtype, 
#         fn_format,
#         sorted,
#         meta_data = [meta_data]
#     )
# end


"""
    defineDataMap(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, String}} = Vector{Dict{String, String}}(),
        constraint::Union{Dict{String, <:AbstractArray{String}}, Nothing} = nothing,
        filename_format::Symbol = :cmip
        dtype::String = "cmip",
        preview::Bool = false,
        sorted::Bool = true
    )

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_data_dirs`, possibly constraint by `constraint`.
"""
function defineDataMap(
    paths_data_dirs::Vector{String}, 
    data_ids::Vector{String};
    constraint::Dict{Symbol, <:AbstractArray{String}} = Dict{Symbol, Vector{String}}(),
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}} = (start_year = typemin(Int), end_year = typemax(Int)),
    level::Symbol = :none,
    dtype::String = "cmip",
    filename_format::Symbol = :cmip,
    sorted::Bool = true,
    #meta_info::Vector{Dict{String, String}} = Dict{String, String}[]
)
    if length(paths_data_dirs) != length(data_ids)
        throw(ArgumentError("'all_paths' and 'ids' must have the same length!"))
    end
    # checkInput with meta_info
    level = toLevel(Val(level))
    fn_format = toFF(Val(filename_format))

    meta_data = _getFilteredMetaData(
        paths_data_dirs, constraint, constraint_ts; level, dtype, fn_format
    )
    _loadDataMapCore(
        meta_data, 
        data_ids;
        constraint_ts,
        level, 
        dtype, 
        fn_format, 
        sorted
        #meta_info
    )
end

function previewDataMap(
    paths_data_dirs::Vector{String}, 
    data_ids::Vector{String};
    constraint::Dict{Symbol, <:AbstractArray{String}} = Dict{Symbol, Vector{String}}(),
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}} = (start_year = typemin(Int), end_year = typemax(Int)),
    level::Symbol = :none,
    dtype::String = "cmip",
    filename_format::Symbol = :cmip
)
    if length(paths_data_dirs) != length(data_ids)
        throw(ArgumentError("'all_paths' and 'ids' must have the same length!"))
    end
    level = toLevel(Val(level))
    fn_format = toFF(Val(filename_format))

    meta_data = _getFilteredMetaData(
        paths_data_dirs, constraint, constraint_ts; level, dtype, fn_format
    )
    if !isa(fn_format, FF_CMIP)
        # if fn type does not have timerange in filename, constraint_ts cannot be applied in preview!
        @info "timeseries constraint not applied for preview when filename format is not :cmip!"
    end
    _previewDataMapCore(meta_data, data_ids)
end


"""
    defineDataMap(
        paths_data_dirs::Vector{String}, 
        data_ids::Vector{String};
        meta_data::Vector{Dict{String, String}} = Vector{Dict{String, String}}(),
        constraint::Union{Dict{String, <:AbstractArray{String}}, Nothing} = nothing,
        filename_format::Symbol = :cmip
        dtype::String = "cmip",
        preview::Bool = false,
        sorted::Bool = true
    )

Return DataMap with entries `data_ids` with the data (model and/or obs depending on `dtype`) 
from all .nc files in all directories in `paths_dirs`, possibly constraint by `constraint`.

For every loaded dataset (entry in built DataMap), files are loaded from several directories;
each entry in `paths_data_dirs` points to the vector of data directories from where data 
is loaded for that dataset.
"""
function defineDataMap(
    paths_data_dirs::Vector{Vector{String}}, 
    data_ids::Vector{String};
    constraint::Dict{String, <:AbstractArray{String}} = Dict{String, Vector{String}}(),
    constraint_ts::Dict{String, Int} = Dict{String, Int}(),
    level::Symbol = :none,
    dtype::String = "cmip",
    filename_format::Symbol = :cmip,
    sorted::Bool = true,
    meta_data::Vector{Dict{String, String}} = Vector{Dict{String, String}}()
)
    level = toLevel(Val(level))
    fn_format = toFF(Val(filename_format))
    checkInput(paths_data_dirs, data_ids; meta_data)
    
    paths_to_files = map(paths_data_dirs) do paths 
        vcat(collectNCFilePaths.(paths)...)
    end
    # if !isnothing(level)        
    #     paths_to_files = filterPathsSharedModels(paths_to_files, level, filename_format)
    # end
    paths_to_files = _applyConstraint(paths_to_files, constraint; level, dtype, fn_format)

    return _loadDataMapCore(
        paths_to_files, 
        data_ids;
        constraint,
        constraint_ts,
        level,
        dtype, 
        fn_format, 
        sorted, 
        meta_data
    )
end


function _getFilteredMetaData(
    paths_data_dirs::Vector{String},
    constraint::Dict{Symbol, <:AbstractArray{String}},
    constraint_ts::NamedTuple{(:start_year, :end_year), <:Tuple{Integer, Integer}};
    level::AbstractLevel = NoLevel(),
    dtype::String = "cmip",
    fn_format::AbstractFnFormat = FF_CMIP()
)
    paths_to_files = collectNCFilePaths.(paths_data_dirs)    
    # TODO: add possibility to have one constraint per dataset
    # first get all the metadata
    meta_data_per_dataset = Vector{Vector{FilenameMeta}}(undef, length(paths_to_files))
    for (i, paths) in enumerate(paths_to_files)
        paths = paths_to_files[i]
        meta = Vector{FilenameMeta}(undef, length(paths))
        for j in eachindex(paths)
            meta[j] = parsePath(paths[j], fn_format)
        end
        meta_data_per_dataset[i] = meta
    end
    # then apply the constraints (TODO: might be done all in once, instead of 2 steps)
    if !isempty(constraint) && any(v -> !isempty(v), values(constraint))
        constraint = Constraint(constraint)
        _applyConstraint!(meta_data_per_dataset, constraint)
    end
    if constraint_ts.start_year != typemin(Int) || constraint_ts.end_year != typemax(Int)
        if isa(fn_format, FF_CMIP)
            # if fn type has timerange in filename, constraint_ts can already be applied!
            _applyConstraintTS!(meta_data_per_dataset, constraint_ts)
        end
    end
    if !isa(level, NoLevel) && dtype == "cmip" # (level doesnt apply to observational data)
        _applyLevel!(meta_data_per_dataset, level)
    end
    return meta_data_per_dataset
end
