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

lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    preview = false # default: false; if true meta data for data to be loaded is returned
);
````
The loaded data is a Vector containing instances of type `Data`. 

We preprocessed the data with ESMValTool using several recipes so that we get separate directories for (the preprocessed data) for every experiment. Thus, to load data for, say lgm and historical experiments, we would call loadDataFromESMValToolConfigs twice with the respective data- and config paths as arguments.

Not all information in the ESMValTool recipes is required to specify the data to be loaded here. For a minimal example of the required structure, see [this example](https://github.com/awi-esc/SimilarityWeights/blob/main/configs/examples/esmvaltool-recipes/mwe_esmvaltool_config.yml).

###  Configuration with seperate yaml file

For an example of a single yaml configuration file [see this example](https://github.com/awi-esc/SimilarityWeights/blob/main/configs/examples/example-lgm-historical.yml).

The config file requires the following entries: 

- `path_data`: base path to where the preprocessed data is stored. This path will be concatenated with the directories specified for each dataset. It can be omitted if inside datasets (see below), base_dir refers to the entire path.

    ```yaml
    path_data: "/albedo/work/projects/p_forclima/preproc_data_esmvaltool"
    ```

- `timerange_to_alias`: mapping from timerange to alias:

    ````yaml
    timerange_to_alias:
    "1850-1900": "historical0"
    "1951-1980": "historical1"
    ````

- `datasets`: list of each dataset to be loaded: 

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

For each given dataset the respective data is loaded from the `base_dir` at
`path_data`. The keys `exp`, `variables` and `statistics` are required to 
specify the data to be loaded. 

`timeranges` is used to filter the data if not 
all available data shall be loaded, but only the data for specific time periods.
`subdirs` is optional, too. If given, data will be loaded only from directories
with names that contain any of the provided values. 



### Subset loaded data

Independently of the way how the data is specified, with a seperate yaml config 
file or using ESMValTool recipes, there's the option to filter the data when 
loading it. Otherwise all  data with the provided properties that is found will
be loaded.

Note that if you specified the data in a new yaml file and further provide 
constraints when loading the data, using the
function `loadDataFromYAML`, the values from the provided function argument 
to filter the data take precedence with respect to what had been specified in 
the yaml file. 


```julia
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes;
    subset = mw.Constraint(
        variables = ["tas", "tos"],
        statistics = ["CLIM"],
        aliases = ["lgm"],
        # timeranges = #  default: ["full"],
        models = Vector{String}(),
        projects = ["CMIP5", "CMIP6"],
        subdirs = ["20241114"]
    )
);
````

The functions `loadDataFromYAML` and `loadDataFromESMValToolConfigs` both 
provide the optional argument `subset` which expects an instance of Type 
`Constraint` which defines the following fields:

- `variables`: the short name of the climate variable, e.g. ["tas", "tos"].

- `statistics`: the statistic/diagnostic that was used when preprocessing the data. The names are generally arbitrary, but need to be identical to those you used in the preprocessing; for instance, we refer to the climatological average as "CLIM", s.t. an example value for this argument is ["CLIM"].

- `aliases`: like for `statistics`, the names are arbitrary but must be identical to what you used when preprocessing the data; e.g. ["historical", "historical1"].

    An alias should refer to a certain `timerange` of a certain experiment, e.g. we call the time period from 1951-1980 'historical1'. To load only this data, it thus does not matter whether you set `timerange=["1951-1980"]` or `alias=["historical1"]` within the `Constraint`-object.

- `timeranges`: timerange; especially for historical data, you may preprocess data for different timeranges,  e.g. ["full", "1980-2014"]. 

- `projects`:  e.g. ["CMIP5", "CMIP6"]. All filenames of the data to be loaded must contain at least one of the given strings. If not specified and `is_model_data=true`, it is set to ["CMIP"]. If not specified and `is_model_data=false`, it is set to ["ERA5"], which is thus the default observational dataset used.

- `models`: can refer to models or individual model members, e.g. ["AWI-CM-1-1-MR"]. All filenames must contain at least one of the given strings + "_". The underscore is important since some models have names that are substrings of other models, e.g. "CNRM-CM5" and "CNRM-CM5-C2". 

- `subdirs`: used to constrain subdirectories to be considered when the data is stored in different subdirectories for each variable (i.e. when `dir_per_var=true`), e.g. ["20241121"] for the case that you only want data from that date and the names of the subdirectories that contain the data include the date.


## How to compute weights

For a full example, see `scripts/run-climwip-simplified.yml`.

Calling `mw.computeWeights(model_data, obs_data, config)` will compute weights based on the provided data.
The config parameter is of type `ConfigWeights` defined in `src/data-utils.jl`. It holds information about the contribution of each combination of statistic/diagnostic and climate variable, once for computing
independence weights and once for computing performance weights. Further parameters from the weighting approach are specified here (`sigmaD`, `sigmaS`). 

The output of the function `computeWeights` is an object of type `ModelWeights` (see `src/data-utils.jl`) which
holds the performance weights for all combinations of statistics/diagnostics and climate variables, for performance as well as independence weights. 
Further, it contains the overall weights (one for each model, summing up to 1) as well as the performance and independence weights 
for each variable (summed across statistics/diagnostics).
<!-- - `weights_variables:`: For each of 'performance' and 'independence' one value per climate variable considered. These values represent the weight of how much each climate variable influences the generalized distance of a model, which is computed by taking a weighted average across the distances with respect to different variables. Should sum up to 1.  -->

