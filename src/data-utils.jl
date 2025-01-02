import YAML
using DimensionalData
using Interpolations


@kwdef struct MetaAttrib
    variable::String
    statistic::String
    alias::String
    exp::String
    timerange::String
end
# Overload the Base.show method to print key-value pairs of MetaAttrib instances
function Base.show(io::IO, x::MetaAttrib)
    println(io, "::$(typeof(x)):")
    for field in fieldnames(MetaAttrib)
        value = getfield(x, field)
        if !isempty(value)
            print(io, "$field=$value ")
        end
    end
end


@kwdef struct Constraint
    variables::Vector{String} = Vector{String}()
    statistics::Vector{String} = Vector{String}()
    aliases::Vector{String} = Vector{String}()
    timeranges::Vector{String} = Vector{String}()
    models::Vector{String} = Vector{String}()
    projects::Vector{String} = Vector{String}()
    subdirs::Vector{String} = Vector{String}()
end


@kwdef struct MetaData
    id::String
    attrib::MetaAttrib
    paths::Vector{String}
end
# Pretty print MetaData instances
function Base.show(io::IO, x::MetaData)
    println(io, "::$(typeof(x))")
    println(io, "$(x.id) (timerange: $(x.attrib.timerange), experiment: $(x.attrib.exp))")
    for path in x.paths
        println(io, "\t$path")
    end
end


@kwdef struct Data
    meta::MetaData
    data::DimArray
end
# Pretty print Data instances
function Base.show(io::IO, x::Data)
    println(io, "::$(typeof(x)) with:")
    print(io, x.meta)
end

@kwdef struct ConfigWeights
    performance::Dict{String, Number}=Dict()
    independence::Dict{String, Number}=Dict()
    sigma_performance::Number=0.5
    sigma_independence::Number=0.5
    ref_period::String=""
    target_path::String=""
end


@kwdef struct Weights
    performance_distances::DimArray
    independence_distances::DimArray
    Di::DimArray # generalized distances each model wrt performance
    Sij::DimArray # generalized distances between pairs of models
    wP::DimArray # normalized
    wI::DimArray # normalized
    w::DimArray # normalized
    w_members::DimArray # weights distributed evenly across resp. model members
    config::ConfigWeights # metadata
end

# Pretty print Weights
function Base.show(io::IO, x::Weights)
    println(io, "Performance Distances ():")
    for m in dims(x.w, :model)
        println(io, "$m: $(round(x.w[model = At(m)], digits=3))")
    end
end

"""
    buildCMIP5EnsembleMember(
        realizations::Vector, initializations::Vector, physics::Vector
    )

Concatenate model settings to build ripf-abbreviations for CMIP5 models which 
do not have it in their metadata.

# Arguments
- `realizations`:
- `initializations`:
- `physics`:
"""
function buildCMIP5EnsembleMember(
    realizations::Vector, initializations::Vector, physics::Vector
)
    n = length(realizations)
    if (n != length(initializations)) || (n != length(physics))
        msg = "inconsistent input for building up CMIP5EnsembleMembers!"
        throw(ArgumentError(msg))
    end
    function concat(elem, prefix)
        if !(elem isa Vector)
            return [prefix * string(elem) for _ in range(1, n)]
        else
            return map(x -> prefix * string(x), elem)
        end
    end
    rips = [concat(realizations, "r"), concat(initializations, "i"), concat(physics, "p")];
    variants = [join(k, "") for k in zip(rips...)]
    return variants
end


"""
    updateMetadata!(
        meta::Dict{String, Union{Array, String}},
        source_names::Vector{String},
        isModelData::Bool
    )

Update metadata 'meta' s.t. attributes that were only present in some files/models 
are set to missing. Further key-value pairs are added concerning the data sources:
For model data:
    - 'member_names': vector that contains for every model a vector with the
    unique names of that model's members
    identifier consisting of variant_label, model_name and for CMIP6 models also grid_label.
    - 'model_names': vector whose length is the sum of the number of all models'
    members; it contains the model names for each unique model member, i.e. this
    vector will not unique if any model had several members
For observational data:
    - 'source_names': vector of data sources

Arguments:
- `meta`:
- `source_names`:
- `isModelData`:
"""
function updateMetadata!(
    meta::Dict{String, Any},
    source_names::Vector{String},
    is_model_data::Bool
)
    sort_indices = sortperm(source_names)
    for key in keys(meta)
        values = meta[key][sort_indices]
        meta[key] = values
        # if none was missing and all have the same value, just use a string
        if !any(ismissing, values) && length(unique(values)) == 1
            meta[key] = string(values[1])
        end
    end
    included_data = Array{String}(source_names[sort_indices])
    if is_model_data
        meta["member_names"] = getUniqueMemberIds(meta, included_data)
        meta["model_names"] = included_data
    else
        meta["source_names"] = included_data
    end
    return nothing
end


"""
    joinMetadata(
    meta1::Dict{String, Any}, meta2::Dict{String, Any}, isModelData::Bool
)

Join the metadata of data to be joint (e.g. when loading data for tos for
CMIP5 and CMIP6 from different locations). If keys are present in 'meta1' or 
'meta2' but not the other, missing values are added.

# Arguments:
- `meta1`:
- `meta2`:
- `is_model_data`: true for model data, false for observational data
"""
function joinMetadata(
    meta1::Dict{String, Any}, meta2::Dict{String, Any}, is_model_data::Bool
)
    meta = Dict{String, Any}()
    n1 = is_model_data ? length(vcat(meta1["member_names"]...)) : length(meta1["source_names"])
    n2 = is_model_data ? length(vcat(meta2["member_names"]...)) : length(meta2["source_names"])
    keys_meta1 = keys(meta1)
    keys_meta2 = keys(meta2)
    keys_shared = collect(intersect(keys_meta1, keys_meta2))
    keys_uniq_m1 = filter(x -> !(x in keys_shared), keys_meta1)
    keys_uniq_m2 = filter(x -> !(x in keys_shared), keys_meta2)

    for k in keys_shared
        v1 = meta1[k]
        v2 = meta2[k]
        if isa(v1, String)
            if isa(v2, String)
                if v1 == v2
                    meta[k] = v1
                else
                    meta[k] = vcat(repeat([v1], outer=n1), repeat([v2], outer=n2))
                end
            else
                meta[k] = vcat(repeat([v1], outer=n1), v2)
            end
        elseif isa(v1, Vector)
            if isa(v2, String)
                meta[k] = vcat(v1, repeat([v2], outer=n2))
            elseif isa(v2, Vector)
                meta[k] = vcat(v1, v2)
            end
        end
    end

    for keys_uniq in [keys_uniq_m1, keys_uniq_m2]
        for k in keys_uniq
            v = get(meta1, k, nothing)
            in_meta1 = true
            if isnothing(v)
                in_meta1 = false
                v = meta2[k]
                n_added = n1
                n = n2
            else
                n_added = n2
                n = n1
            end
            v_added = repeat([missing], outer=n_added)
            # be sure to add vectors in correct order (because v may refer to value of meta1 or meta2!)
            if isa(v, String)
                if in_meta1
                    meta[k] = vcat(repeat([v], outer=n), v_added)
                else
                    meta[k] = vcat(v_added, repeat([v], outer=n))
                end
            else
                if in_meta1
                    meta[k] = vcat(v, v_added)
                else
                    meta[k] = vcat(v_added, v)
                end
            end
        end
    end
    return meta
end


"""
    updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)

Vectors in metadata 'meta' refer to different models (members).
These are now summarized such that each vector only contains N entries where N
is the number of models (i.e. without the unique members).
If the metadata for members of a model differ across members, the respective
entry in the vector will be a vector itself.

# Arguments:
- `meta`:
- `grouped_data`:
"""
function updateGroupedDataMetadata(meta::Dict, grouped_data::DimensionalData.DimGroupByArray)
    meta_new = filter(((k,v),) -> k=="member_names" || !(v isa Vector), meta)
    attributes = filter(x -> meta[x] isa Vector && x != "member_names", keys(meta))
    attribs_diff_across_members = [];
    # iterate over attributes that are vectors, thus different for the different
    # members or models
    for key in attributes
        for (i, model) in enumerate(dims(grouped_data, :model))
            indices = findall(x -> x==model, meta["model_names"])
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
        @warn "metadata attributes that differ across model members (ok for some!)" unique(attribs_diff_across_members)
    end
    return meta_new
end


"""
    computeInterpolatedWeightedQuantiles(
        quantiles::Vector{<:Number},
        vals::Vector;
        weights=nothing    
    )

This implementation follows the one used by Brunner et al.

# Arguments:
- `quantiles`: 
- `vals`:
- `weights`:
"""
function computeInterpolatedWeightedQuantiles(
    quantiles::Vector{<:Number},
    vals::Vector;
    weights=nothing
)
    if isnothing(weights)
        weights = ones(length(vals));
    end
    indicesSorted = Array(sortperm(vals)); # gives indices of smallest to largest data point
    weightsSorted = weights[indicesSorted];
    weightedQuantiles = cumsum(weightsSorted) - 0.5 * weightsSorted;
    weightedQuantiles = reshape(weightedQuantiles, length(weightedQuantiles), 1);
    weightedQuantiles = (weightedQuantiles .- minimum(weightedQuantiles)) ./ maximum(weightedQuantiles);

    interp_linear = Interpolations.linear_interpolation(
        vec(weightedQuantiles),
        vals[indicesSorted],
        extrapolation_bc=Interpolations.Line()
    );
    return interp_linear(quantiles)
end


"""
    setLookupsFromMemberToModel(data::DimArray, dim_names::Vector{String})

Change the lookup values for the dimension 'member' to refer to the models, i.e.
they are not unique anymore. This is done in preparation to group the data by
the different models.

# Arguments:
- `data`: has at least dimensions in 'dim_names'
- `dim_names`: names of dimensions to be changed, e.g. 'member', 'member1'
"""
function setLookupsFromMemberToModel(data::DimArray, dim_names::Vector{String})
    n_dims = length(dim_names)
    for (i, dim) in enumerate(dim_names)
        unique_members = dims(data, Symbol(dim))
        models = map(x -> split(x, MODEL_MEMBER_DELIM)[1], unique_members)

        data = DimensionalData.set(data, Symbol(dim) => models)
        new_dim_name = n_dims > 1 ? "model" * string(i) : "model"
        data = DimensionalData.set(data, Symbol(dim) => Symbol(new_dim_name))
    end
    return data
end


"""
    getMetaAttributesFromESMValToolConfigs(
        base_path_configs::String;
        subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
)

Read config files (ESMValTool recipes) to determine which data shall be loaded
in general (datapaths are not handled here), as defined by variable, statistic, 
experiment and timerange/alias.

# Arguments:
- `base_path_configs`: paths to directory that contain one or more yaml configs,
which may be ESMValTool recipes.
- `subset`:

# Returns: A Vector of `MetaAttrib`-Objects.
"""
function getMetaAttributesFromESMValToolConfigs(
    base_path_configs::String;
    constraint::Union{Constraint, Nothing} = nothing
)
    paths_configs = filter(
        x -> isfile(x) && endswith(x, ".yml"),
        readdir(base_path_configs, join=true)
    )
    meta_attribs = Vector{MetaAttrib}()
    for path_config in paths_configs
        config = YAML.load_file(path_config)
        data_all = config["diagnostics"]
        aliases = keys(data_all)

        for alias in aliases
            data = data_all[alias]["variables"]
            for (k,v) in data
                variable, statistic = String.(split(k, "_"))
                if typeof(v["exp"]) <: String
                    experiment = v["exp"]
                else
                    experiment = join(v["exp"], "-")
                end
                timerange = replace(get(v, "timerange", "full"), "/" => "-")
                meta = MetaAttrib(
                    variable = variable, 
                    statistic = statistic, 
                    alias = alias, 
                    exp = experiment,
                    timerange = timerange
                )
                if !(meta in meta_attribs) 
                    push!(meta_attribs, meta)
                end
            end
        end
    end
    if !isnothing(constraint)
        applyDataConstraints!(meta_attribs, constraint)
    end
    return meta_attribs
end


function addMetaData!(meta_dict::Dict{String, MetaData}, meta::MetaData)
    id = getMetaDataID(meta.attrib)
    if haskey(meta_dict, id)
        paths = meta_dict[id].paths
        append!(paths, meta.paths)
        meta_dict[id] = MetaData(
            id = id,
            attrib = meta.attrib,
            paths = unique(paths)
        )
    else
        meta_dict[id] = meta
    end
    return nothing
end

"""
    getMetaDataFromYAML(
    path_config::String,
    dir_per_var::Bool,
    is_model_data::Bool
    subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
)

Load data as specified in config file located at `path_config`. For constraints
that are specified in the config file as well as in the `subset` argument, 
the values of the latter have precedence over the former. The constraints given
in the argument `subset` are applied to EVERY dataset specified in the config 
file.

# Arguments:
- `path_config`: path to config yaml file specifying meta attributes and pathes of data
- `dir_per_var`: true if data for each climate variable is stored in seperate directory 
- `is_model_data`: true for model data, false for observational data
- `subset`: TODO
"""
function getMetaDataFromYAML(
    path_config::String, 
    dir_per_var::Bool, 
    is_model_data::Bool;
    constraint::Union{Constraint, Nothing} = nothing
)
    config = YAML.load_file(path_config)
    datasets = config["datasets"]
    base_path = get(config, "path_data", "")
    timerange_to_alias = config["timerange_to_alias"]

    """
    If data is constraint by provided argument when loading, the argument takes
    precedence over the given value inside the config yaml file.
    """
    function getConstraintVal(
        cs_arg::Union{Constraint, Nothing}, 
        cs_config::Dict, 
        field::String;
        default_val::Vector{String} = Vector{String}()
    )
        result = get(cs_config, field, default_val)
        if !isnothing(cs_arg)
            value = getfield(cs_arg, Symbol(field))
            result = isempty(value) ? result : value
        end
        return result
    end

    meta_data = Dict{String, MetaData}()
    for ds in datasets
        data_dir = get(ds, "base_dir", nothing)
        experiment = get(ds, "exp", nothing)
        if isnothing(experiment) || isnothing(data_dir)
            msg = "Config yaml file must specify values for keys 'exp' (experiment) and 'base_dir' (path to data directory)!"
            throw(ArgumentError(msg))
        end
        variables = getConstraintVal(constraint, ds, "variables")
        statistics = getConstraintVal(constraint, ds, "statistics")
        if isnothing(variables) || isnothing(statistics)
            throw(ArgumentError("Config yaml file must specify values for keys 'variables', 'statistics'!"))
        end
        path_data = joinpath(base_path, data_dir)
        # timeranges is optional in config file; they can be subset by aliases provided in constraint argument
        # in config file only timeranges are provided/considered, not aliases
        timeranges = getConstraintVal(constraint, ds, "timeranges"; default_val = ["full"])
        aliases = isnothing(constraint) ? Vector{String}() : constraint.aliases
        if !isempty(aliases)
            # use subset of aliases as given by input argument constraint
            filter!(x -> get(timerange_to_alias, x, nothing) in constraint.aliases, timeranges)
            if isempty(timeranges)
                msg = "Given aliases in constraint do not match with considered timeranges for dataset: $ds !"
                throw(ArgumentError(msg))
            end
        end
        projects = getConstraintVal(constraint, ds, "projects")
        models = getConstraintVal(constraint, ds, "models")
        subdirs = getConstraintVal(constraint, ds, "subdirs")
        cs = Constraint(
            statistics=statistics, 
            aliases=aliases, 
            timeranges=timeranges,
            models=models, 
            projects=projects, 
            subdirs=subdirs
        )

        for clim_var in variables
            for stat in statistics
                for timerange in timeranges
                    alias = timerange == "full" ? experiment : 
                        get(timerange_to_alias, timerange, "unknown")
                    attribs = [MetaAttrib(
                        variable = clim_var, 
                        statistic = stat, 
                        alias = alias, 
                        exp = experiment,
                        timerange = timerange
                    )]
                    # apply constraints if given
                    constraint_vals = map(x->getfield(cs, x), fieldnames(Constraint))
                    if any(!isempty, constraint_vals)
                        applyDataConstraints!(attribs, cs)
                    end
                    if !isempty(attribs)
                        meta = buildMetaData(
                            attribs[1], path_data, dir_per_var, is_model_data; 
                            constraint=cs
                        )
                        addMetaData!(meta_data, meta)
                    end
                end
            end
        end
    end
    return meta_data
end


"""
    buildPathsToDataFiles(
        path_data::String,
        is_model_data::Bool;
        subset::Union{Dict{String, Vector{String}}, Nothing}=nothing
    )

# Arguments:
- `path_data`:
- `is_model_data`: true for model data false for observational data
- `subset`:

# Returns: Vector of Strings containing paths to data files in `path_data` 
that were not filtered out by `subset`.
"""
function buildPathsToDataFiles(
    path_data::String,
    is_model_data::Bool;
    model_constraints::Vector{String} = Vector{String}(),
    project_constraints::Vector{String} = Vector{String}()
)
    if !isdir(path_data)
        throw(ArgumentError(path_data * " does not exist!"))
    end
    if isempty(project_constraints)
        project_constraints = is_model_data ? ["CMIP"] : ["ERA5"]
    end
    ncFiles = filter(
        x -> isfile(x) && endswith(x, ".nc"),
        readdir(path_data; join=true)
    )
    paths_to_files = Vector{String}()
    # constrain files that will be loaded
    for file in ncFiles
        keep = !isempty(project_constraints) ? 
            any([occursin(name, file) for name in project_constraints]) : true
        if !keep
            @debug "exclude $file because of projects subset"
            continue
        end
        keep = !isempty(model_constraints) ? 
            applyModelConstraints(file, model_constraints) : true
        if keep
            push!(paths_to_files, file)
        end
    end
    return paths_to_files
end


function getMetaDataID(attrib::MetaAttrib)
    # id=join([getfield(attrib, field) for field in fieldnames(MetaAttrib)], "_")
    # return id
    return join([attrib.variable, attrib.statistic, attrib.alias], "_")
end


"""
    buildMetaData(
        attrib::Union{MetaAttrib, Vector{MetaAttrib}},
        base_path_data::String,
        dir_per_var::Bool,
        is_model_data::Bool;
        subset::Union{Dict{String, Vector{String}}, Nothing} = nothing
    )

Create data paths from assumed underlying data structure and `base_path_data` 
for every `MetaAttrib`-object in `attrib` and return a Vector of `MetaData`-
objects.

# Arguments:
- `attrib`:
- `base_path_data`:
- `dir_per_var`:
- `is_model_data`:
- `subset`:

# Return: Vector of `MetaData`-objects.
"""
function buildMetaData(
    attrib::MetaAttrib,
    base_path_data::String,
    dir_per_var::Bool,
    is_model_data::Bool;
    constraint::Union{Constraint, Nothing} = nothing
)
    subdir_constraints = isnothing(constraint) ? nothing : constraint.subdirs
    paths_data = buildPathsForMetaAttrib(
        base_path_data, attrib, dir_per_var; subdir_constraints
    )
    paths_to_files = Vector{String}()
    for path_data in paths_data
        paths = buildPathsToDataFiles(
            path_data, is_model_data; 
            model_constraints=constraint.models,
            project_constraints = constraint.projects
        )
        append!(paths_to_files, paths)
    end
    # for observational data, experiment doesn't make sense
    if !is_model_data
        attrib = MetaAttrib(
            variable = attrib.variable, 
            statistic = attrib.statistic,
            alias = attrib.alias,
            exp = "",
            timerange = attrib.timerange
        )
    end
    return MetaData(
        id = getMetaDataID(attrib),
        attrib = attrib,
        paths = paths_to_files
    )    
end

# if dir_per_var is true, directories at base_paths have subdirectories,
# one for each variable (they must contain '_VAR', e.g. '_tas'),
# otherwise base_paths are the paths to the directories that contain
# a subdirectory 'preproc'
"""
    buildPathsForMetaAttrib(
        base_path::String, 
        attrib::MetaAttrib,
        dir_per_var::Bool;
        subdir_constraints::Union{Vector{String}, Nothing}=nothing
    )

# Arguments:
- `base_path`: base directory of stored data specified in `attrib`.
- `attrib`: meta attributes of data.
- `dir_per_var`: true if data of each climate variable is stored in a seperate directory.
- `subdir_constraints`: if given, paths must contain ANY of the given elements. Existing paths that don't are ignored.
"""
function buildPathsForMetaAttrib(
    base_path::String, 
    attrib::MetaAttrib,
    dir_per_var::Bool;
    subdir_constraints::Union{Vector{String}, Nothing}=nothing
)
    base_paths = [base_path]
    if dir_per_var
        base_paths = filter(isdir, readdir(base_path, join=true))
        filter!(x -> occursin("_" * attrib.variable, x), base_paths)
    end
    if !isnothing(subdir_constraints) && !isempty(subdir_constraints)
        filter!(p -> any([occursin(name, p) for name in subdir_constraints]), base_paths)
    end
    data_paths = Vector{String}()
    for path in base_paths
        # Note: particular data structure assumed here!
        diagnostic = join([attrib.variable, attrib.statistic], "_")
        path_data = joinpath(path, "preproc", attrib.alias, diagnostic)
        if !isdir(path_data)
            @warn "No data found at:" path_data
        else
            push!(data_paths, path_data)
        end
    end
    return data_paths
end


"""
    applyDataConstraints!(
    meta_attributes::Vector{MetaAttrib}, subset::Dict{String, Vector{String}}
)

Subset model ids so that only those with properties specified in 'subset' remain.

# Arguments
- `meta_attributes`: Vector of `MetaAttrib` instances.
- `subset`: Mapping from fieldnames of 'MetaAttrib' struct to Vector 
specifiying the properties of which at least one must be present for an id to 
be retained.
"""
function applyDataConstraints!(
    meta_attributes::Vector{MetaAttrib}, constraint::Constraint
)
    # check for compatibility of timerange and alias first
    timerange_constraints = getfield(constraint, :timeranges)
    alias_constraints = getfield(constraint, :aliases)
    timerangeOk(attrib::MetaAttrib) = any(x -> attrib.timerange == x, timerange_constraints)
    aliasOk(attrib::MetaAttrib) = any(x -> attrib.alias == x, alias_constraints)

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

    # constraints wrt variables and statistics
    stats_constraints = constraint.statistics
    vars_constraints = constraint.variables
    stats_ok(attrib::MetaAttrib) = any(x -> attrib.statistic == x, stats_constraints)
    vars_ok(attrib::MetaAttrib) = any(x -> attrib.variable == x, vars_constraints)
    if !isempty(stats_constraints)
        filter!(stats_ok, meta_attributes)
    end
    if !isempty(vars_constraints)
        filter!(vars_ok, meta_attributes)
    end
    return nothing
end


function applyModelConstraints(file::String, model_constraints::Vector{String})
    # model constraints may contain individual members
    # (e.g. for "CNRM-CM5#r1i1p1", the model name, CNRM-CM5, as well as the
    # member id, r1i1p1, have to be part of the filename,
    # but not with the delimiter # as given here)
    # for CMIP6 models we further added the grid to the member id at the end 
    # after an underscore, i.e. for CMIP6 it is assumed that the filename 
    # ends with _GRID, e.g. _gn.nc
    split_chars = Regex("[$(MODEL_MEMBER_DELIM)_]")
    model_member_constraints = map(x -> split(x, split_chars), model_constraints)
    keep_file = true
    for constraints in model_member_constraints
        # adding the suffix "_" is important since otherwise, for instance,
        # CNRM-CM5-C2 would remain even if the constraint was a substring like
        # CNRM-CM5
        keep_file = all([occursin(name * "_", file) for name in constraints[1:2]])
        if keep_file && length(constraints) == 3 # CMIP6 with grid
            keep_file = occursin("_" * constraints[3], file)
        end
        if keep_file
            break
        end
    end
    if !keep_file
        @debug "exclude $file because of model constraints (or shared models): $model_constraints"
    end
    return keep_file
end


"""
    indexData(data::Data, clim_var::String, diagnostic::String, alias::String)
"""
function indexData(data::Vector{Data}, clim_var::String, diagnostic::String, alias::String)
    fn(data::Data) = 
        data.meta.attrib.variable == clim_var && 
        data.meta.attrib.statistic == diagnostic &&
        data.meta.attrib.alias == alias
    
    df = filter(fn, data)
    if length(df) == 0
        throw(ArgumentError("No data contained for $clim_var, $diagnostic, $alias !"))
    elseif length(df) > 1
        @warn "more than one dataset given for $clima_var, $diagnostic, $alias."
    end
    return df[1].data
end


"""
    computeDistancesAllDiagnostics(
        model_data::Vector{Data}, 
        obs_data::Union{Nothing, Data}, 
        config::Dict{String, Number},
        ref_period::String,
        forPerformance::Bool    
    )

Compute RMSEs between models and observations or between predictions of models.

# Arguments:
- `model_data`:
- `obs_data`:
- `config`:
- `ref_period`:
- `for_performance`: true for distances between models and observations, false for distances between model predictions
"""
function computeDistancesAllDiagnostics(
    model_data::Vector{Data}, 
    obs_data::Union{Nothing, Vector{Data}}, 
    config::Dict{String, Number},
    ref_period_alias::String,
    for_performance::Bool
)
    # compute performance/independence distances for all model members
    var_diagnostic_keys = collect(keys(config))
    diagnostics = String.(unique(map(x -> split(x, "_")[2], var_diagnostic_keys)))
    distances_all = []
    for diagnostic in diagnostics
        distances = []
        diagnostic_keys = filter(x -> endswith(x, "_" * diagnostic), var_diagnostic_keys)
        variables = String.(map(x -> split(x, "_")[1], diagnostic_keys))
        for clim_var in variables
            models = indexData(model_data, clim_var, diagnostic, ref_period_alias)

            if for_performance
                observations = indexData(obs_data, clim_var, diagnostic, ref_period_alias)
                if length(dims(observations, :source)) != 1
                    @warn "several observational datasets available for computing distances"
                end
                observations = observations[source=1]
                dists = getModelDataDist(models, observations)
            else
                dists = getModelDistances(models)
            end
            push!(distances, dists)
        end
        # TODO: recheck the metadata here! standard_name should be converted 
        # to a vector?!
        distances = cat(distances..., dims = Dim{:variable}(String.(collect(variables))));
        push!(distances_all, distances)
    end
    return cat(distances_all..., dims = Dim{:diagnostic}(String.(collect(diagnostics))));
end


"""
    computeGeneralizedDistances(
        distances_all::DimArray, weights::DimArray, for_performance::Bool
)

# Arguments:
- `distances_all`:
- `weights`:
- `for_performance`: 
"""
function computeGeneralizedDistances(
    distances_all::DimArray, weights::DimArray, for_performance::Bool
)
    dimensions = for_performance ? (:member,) : (:member1, :member2)
    norm = mapslices(Statistics.median, distances_all, dims=dimensions)
    normalized_distances =  DimArray(
        distances_all ./ norm, dims(distances_all), metadata = distances_all.metadata
    )
    distances = for_performance ? 
        summarizeEnsembleMembersVector(normalized_distances, false) :
        averageEnsembleMatrix(normalized_distances, false);

    distances = mapslices(x -> x .* weights, distances, dims=(:variable, :diagnostic))
    return dropdims(
        sum(distances, dims=(:variable, :diagnostic)), dims=(:variable, :diagnostic)
    )
end


"""
    allcombinations(v...)

Generate all possible combinations of input vectors, where each combination consists of one element from each input vector, concatenated as a string with underscores separating the elements.

# Arguments
- `v...`: A variable number of input vectors.

# Returns
A vector of strings, where each string represents a unique combination of elements from the input vectors, joined by underscores.

# Example
```jldoctest
julia> ModelWeights.allcombinations(["tos", "tas"], ["CLIM"])
2-element Vector{String}:
 "tos_CLIM"
 "tas_CLIM"
```
"""
function allcombinations(v...)
    combis = Vector{String}()
    for elems in Iterators.product(v...)
        push!(combis, join(elems, "_"))
    end
    return combis
end


"""
    isValidDataAndWeightInput(
        data::Data, keys_weights::Vector{String}, ref_period::String
    )

Check that there is data for all keys of the provided (balance) weights and the given reference period 'ref_period'.

# Arguments:
- `data`:
- `keys_weights`:
- `ref_period_alias`:
"""
function isValidDataAndWeightInput(
    data::Vector{Data}, keys_weights::Vector{String}, ref_period_alias::String
)
    actual_ids_data = map(x -> x.meta.id, data)
    required_keys_data = map(x -> x * "_" * ref_period_alias, keys_weights)
    return all([k in actual_ids_data for k in required_keys_data])
end


"""
    getTimerangeAsAlias(meta_attribs::Vector{MetaAttrib}, timerange::String)

Translate given timerange to corresponding alias in 'data_ids'.

# Arguments:
- `meta_attribs`:
- `timerange`:
"""
function getTimerangeAsAlias(meta_attribs::Vector{MetaAttrib}, timerange::String)
    attribs = filter(x -> x.timerange == timerange, meta_attribs)
    return isempty(attribs) ? nothing : attribs[1].alias
end


"""
    getAliasAsTimerange(meta_attribs::Vector{MetaAttrib}, alias::String)

Translate given alias for every `MetaAttrib` in `meta_attribs`to corresponding timerange.

# Arguments:
- `meta_attribs`:
- `alias`:
"""
function getAliasAsTimerange(meta_attribs::Vector{MetaAttrib}, alias::String)
    attribs = filter(x -> x.alias == alias, meta_attribs)
    return isempty(attribs) ? nothing : attribs[1].timerange
end

"""
    getRefPeriodAsTimerangeAndAlias(
        meta_attribs::Vector{MetaAttrib}, ref_period::String
    )
"""
function getRefPeriodAsTimerangeAndAlias(
    meta_attribs::Vector{MetaAttrib}, ref_period::String
)
    alias = getTimerangeAsAlias(meta_attribs, ref_period) 
    # true : ref_period given as alias, false: ref_period given as timerange
    timerange = isnothing(alias) ? getAliasAsTimerange(meta_attribs, ref_period) : ref_period
    alias = isnothing(alias) ? ref_period : alias

    if (isnothing(alias) || isnothing(timerange))
        throw(ArgumentError("ref period $ref_period not given in data!"))
    end
    return (alias = alias, timerange = timerange)
end

"""
    setRefPeriodInWeightsMetadata!(meta::Dict, alias::String, timerange::String)
"""
function setRefPeriodInWeightsMetadata!(meta::Dict, alias::String, timerange::String)
    meta["ref_period_alias"] = alias
    meta["ref_period_timerange"] = timerange
    return nothing
end