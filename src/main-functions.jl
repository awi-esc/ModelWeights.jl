import YAML
using CairoMakie

@kwdef struct Config 
    base_path::String
    target_dir::String
    experiment::String
    prefix_var_folders::String = ""
    variables::Vector{String}
    name_ref_period::String
    models_project_name::String
    obs_data_name::String
    weights_variables::Dict{String, Dict{String,Number}}
    weight_contributions::Dict{String, Number}
end


function buildPathsToVarData(config::Config)
    prefix = config.prefix_var_folders;
    var_to_path = Dict{String, String}();
    for var in config.variables
        folder = ifelse(isempty(prefix), var, prefix * "_" * var);
        path_data = joinpath(
            config.base_path,
            config.experiment, 
            folder,
            "preproc"
        )
        var_to_path[var] = joinpath(
            path_data, 
            "climatology_" * config.name_ref_period,
            var
        )
    end
    return var_to_path
end



function validateConfig(path_config::String)
    config_yaml = YAML.load_file(path_config);
    config = Config(
        base_path = config_yaml["base_path"],
        target_dir = config_yaml["target_dir"],
        experiment = config_yaml["experiment"],
        prefix_var_folders = config_yaml["prefix_var_folders"],
        variables = config_yaml["variables"],
        name_ref_period = config_yaml["name_ref_period"],
        models_project_name = config_yaml["models_project_name"],
        obs_data_name = config_yaml["obs_data_name"],
        weights_variables = convert(
            Dict{String, Dict{String, Number}}, 
            config_yaml["weights_variables"]
        ), 
        weight_contributions = config_yaml["weight_contributions"]
    )
    # TODO: add checks consistency, paths for all specified variables there and exist, etc.
    return config
end

function buildFullModelNames()
end


"""
    getCommonModelsAcrossVars(modelData::Dict{String, DimArray})

Keep only those models for which there is data for all variables.

# Arguments
- modelData: maps from climate variable to respective data.
"""
function getCommonModelsAcrossVars(modelData::Dict{String, DimArray})
    data_all = deepcopy(modelData);
    variables = keys(modelData);
    
    shared_models =  nothing
    for var in variables
        data_var = modelData[var];
        names = Array(dims(data_var, :model));
        variants = data_var.metadata["variant_label"];
        grids = data_var.metadata["grid_label"];
        
        models = map(
            x->join(x, "__", "__"), 
            zip(names, variants, grids)
        );
        data_all[var].metadata["full_model_names"] = models;
        if isnothing(shared_models)
            shared_models = models
        else
            shared_models = intersect(shared_models, models)
        end
    end
    for var in variables
        data = data_all[var]
        indices = findall(m -> m in shared_models, data.metadata["full_model_names"]);
        data_all[var] = data[model = indices];
        # adapt metadata accordingly
        data.metadata["full_model_names"] = data.metadata["full_model_names"][indices];
        for key in SimilarityWeights.LOCAL_METADATA_KEYS
            data.metadata[key] = data.metadata[key][indices];
        end
    end
    # sizes = map(v-> size(data_all[v]), collect(variables));
    # if length(unique(sizes)) != 1
    #     @warn "Models were not uniquely identified by "
    # end
    return data_all
end


"""
    run(path_config::String, plot::Bool=false)

Note that the combined weights can have different dimension of models since 
all ensemble members of th same model are averaged before combining performance
and independence weights.
"""
function runWeights(path_config::String, plot::Bool=false)
    config = validateConfig(path_config);
    pathsDict = buildPathsToVarData(config)
    modelData = loadPreprocData(pathsDict, [config.models_project_name]);
    obsData = loadPreprocData(pathsDict, [config.obs_data_name])
    
    modelDataAllVars =  getCommonModelsAcrossVars(modelData);
    wP = getPerformanceWeights(modelDataAllVars, obsData, config.weights_variables["performance"]);
    wI = getIndependenceWeights(modelDataAllVars, config.weights_variables["independence"]);
    sigmas = config.weight_contributions
    weights = combineWeights(wP, wI, sigmas["performance"], sigmas["independence"]);

    result = Dict(
        "performance" => wP,
        "independence" => wI,
        "combined" => weights
    );

    if plot
        for var in keys(modelDataAllVars)
            means = dropdims(mean(modelDataAllVars[var], dims=:model), dims=:model);
            title = join(["Unweighted average;", var, "in", 
                means.metadata["units"], "\n experiment:", means.metadata["experiment_id"]], " ", " ");
            target = Target(
                directory = config.target_dir,
                filename = "unweighted_avg_" * var * ".png",
                save = true
            )
            fig = SimilarityWeights.plotMeansOnMap(means,  title, target);

        end
    end
    return result
end

