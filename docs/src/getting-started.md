# Getting started

In the following, we explain how to load data and how to compute model weights.

## How to load data

The data to be loaded is specified in yaml config files. There are two 
possibilities how to organize these config files: either we use one or more 
ESMValTool recipes that we had used to preprocess the data or we write a new 
yaml file, independently of the ESMValTool recipes.

### Configuration with ESMValTool recipes

We use the completely spelled out versions of the recipes that ESMValTool
returns (RECIPENAME_filled.yml, stored in run-directory).

````julia
import ModelWeights as mw

# set path_data to directory where your preprocessed data is stored
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
# set path_recipes: directory where the ESMValTool recipes are stored 
path_recipes = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/configs-ModelWeights/esmvaltool-recipes/lgm-cmip5-cmip6";

lgm_data = mw.loadDataFromESMValToolRecipes(
    path_data, path_recipes;
    dir_per_var = true, # default: true
    is_model_data = true, # default: true
    subset = nothing, # default: nothing
    preview = false # default: false; if true meta data for data to be loaded is returned
);
````
The loaded data is a Dictionary mapping from an identifier of the form 'variable\_diagnostic\_alias' (e.g., tas\_CLIM\_lgm) to a `YAXArray`. Additional metadata that was added by us is part of the YAXArrays .properties dictionary. All such fields start with an underscore, e.g. '\_id', '\_paths', '\_statistic'.

Not all information in the ESMValTool recipes is required to specify the data to be loaded here. For a minimal example of the required structure, see [this example](https://github.com/awi-esc/ModelWeights/blob/main/configs/examples/esmvaltool-recipes/mwe_esmvaltool_config.yml).

```@raw html
<!-- TODO: explain dir_per_var argument -->
```
How the loaded data can be constrained by the subset argument of the function `loadDataFromESMValToolRecipes` is described further below. 

###  Configuration with seperate yaml file
For an example of a single yaml configuration file [see this example](https://github.com/awi-esc/ModelWeights/blob/main/configs/examples/example-lgm-historical.yml).

Keys expected in the config-yaml file:

- `path_data`: base path to where the preprocessed data is stored. This path will be concatenated with the directories specified for each dataset. It can 
be omitted if inside datasets (see below), `base_dir` refers to the entire absolute path.

- `datasets`: list of each dataset to be loaded. Each dataset is specified by a dictionary with the following keys:
    
    - `base_dir::String`: relative path to the data with respect to the path given in `path_data`. If `path_data` is not specified in the config file, `path_data` for this dataset is assumed to refer to an absolute path instead.
    - `exp::String`: e.g., "historical"
    - `variables::Vector{String}`: e.g., ["tos", "tas"]
    - `statistics::Vector{String}`: may be empty for fixed variables for which no statistics are computed. In that case, `statistics` is set to ["none"] and a warning is thrown unless for known fixed variables (for now just [`orog`]).

   The following keys are optional:
   - `subdirs`: vector of strings; only subdirectories that contain any of the given strings in their name will be considered.
   - `timeranges`:
   - `aliases`:

````yaml
datasets: [
{
    base_dir: "LGM", 
    exp: "lgm", 
    variables: ["tas", "tos"], 
    statistics: ["CLIM"], 
    subdirs: ["20241114"]
},
{
    base_dir: "historical", # required String
    exp: "historical", # required String
    variables: ["tas", "tos"], # required Vector
    statistics: ["CLIM"], # required Vector
    timeranges: ["full"], # optional Vector
    subdirs: ["20241121", "20241118"] # optional Vector
}
]
````

- `timerange_to_alias`: mapping from timerange to alias.

````yaml
timerange_to_alias:
"1850-1900": "historical0"
"1951-1980": "historical1"
````


For each given dataset the respective data is loaded from the `base_dir` at
`path_data`. The keys `base_dir`, `exp`, `variables` and `statistics` are required to 
specify the data to be loaded. 
`timeranges` and `subdirs` are optional and work in the same way as when given as field of the optional argument `subset` described below.

As we'll explain next, you can also provide further constraints when loading the data using the function `loadDataFromYAML`. Note that in this case, the values from the provided function argument (subset::Constraint) to filter the data take precedence over what had been specified in the yaml file. 


### Optional parameters for filtering data
For both functions, `loadDataFromESMValToolRecipes` and `loadDataFromYAML`, there is a set of optional parameters in order to constrain the loaded data:

- `preview`: If set to false (default), the data will not be loaded and only the metadata with the information that we added that specifies which data will be loaded is returned.

- `dir_per_var`: If set to true (default), only subdirectories of the base\_path that contain `_VARIABLE` in their name will be searched. 

- `is_model_data`: is set to true (default) when loading CMIP data, to false when loading observationa/reanalysis data (e.g. ERA5).

```@raw html
<!-- `subset_shared`: Can be either `nothing` (default), `ModelWeights.MODEL` or `ModelWeights.MEMBER`. If set to `MODEL`, only data from the same models will be loaded. If, for instance, data from lgm and historical experiments shall be loaded, this configuration will ensure that for both, only the same models are loaded, i.e. those for which both experiments were done. While this considers models, not specific simulations, setting subset\_shared to `MEMBER` would only load models that share the exact same simulations (i.e. the same member\_id abbreviation, e.g. `r1i1p1f1`). -->
```

- `subset` is an optional parameter of type `Constraint` or `Nothing`. A Constraint further constrains the data to be loaded and has the following fields:
   
    - `variables`: The short name of the climate variable, e.g. ["tas", "tos"].
    
    - `statistics`: The statistic/diagnostic that was used when preprocessing the data. The names are generally arbitrary, but need to be identical to those you used in the preprocessing; for instance, we refer to the climatological average as "CLIM", s.t. an example value for this argument is ["CLIM"].
    
    - `aliases`: Like for `statistics`, the names are arbitrary but must be identical to what you used when preprocessing the data; e.g. ["historical", "historical1"].
    An alias should refer to a certain `timerange` of a certain experiment, e.g. we call the time period from 1951-1980 'historical1'. To load only this data, it thus does not matter whether you set `timerange=["1951-1980"]` or `alias=["historical1"]`.
    
    - `timeranges`: Timeranges to be loaded; especially for historical data, you may preprocess data for different timeranges,  e.g. ["full", "1980-2014"]. 

    - `projects`:  e.g. ["CMIP5", "CMIP6"]. All filenames of the data to be loaded must contain at least one of the given strings. If not specified and `is_model_data=true`, it is set to ["CMIP"]. If not specified and `is_model_data=false`, it is set to ["ERA5"], the default observational dataset used.

    - `models`: List of models or individual model members, e.g. ["AWI-CM-1-1-MR"]. All filenames must contain at least one of the given strings + "_". The underscore is important since some models have names that are substrings of other models, e.g. "CNRM-CM5" and "CNRM-CM5-C2". 
   
    - `subdirs`: If given, data will be loaded only from subdirectories of the given base_dir that contain any of the provided values in their name. This is recommended when there are many subdirectories for a specific variable  within base\_dir and you only want data from a specific one (e.g. of a certain date, given that the date is included in the name of the directory).



## How to compute weights

For a full example, see `scripts/run-climwip-simplified.yml`.

```@raw html
<!-- TODO: update the following -->
```
Calling `mw.computeWeights(dists_indep, dists_perform, config)` will compute weights based on the distances for the independence weights and the distances for the performance weights.
The config parameter is of type `ConfigWeights` defined in `src/data-utils.jl`. It holds information about the contribution of each combination of statistic/diagnostic and climate variable, once for computing independence weights and once for computing performance weights. Further parameters concerning the weighting are specified in the ConfigWeights struct, such as the hyperparameters,  `sigma_performance` and `sigma_independence`. 

```@raw html
<!-- TODO: update following-->
```
The output of the function `computeWeights` is an object of type `Weights` (see `src/data-utils.jl`) which holds the weights for all combinations of statistics/diagnostics and climate variables, for performance as well as independence weights. Further, it contains the overall weights (one for each model, summing up to 1) as well as the performance and independence weights for each variable (summed across statistics/diagnostics).

```@raw html
<!-- `weights_variables:`: For each of 'performance' and 'independence' one value per climate variable considered. These values represent the weight of how much each climate variable influences the generalized distance of a model, which is computed by taking a weighted average across the distances with respect to different variables. Should sum up to 1.  -->
```
