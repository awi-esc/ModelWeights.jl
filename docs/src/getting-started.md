# Getting started
## How to load data

Data is loaded with the function `loadData` which has two positional arguments: `base_path` which is a String pointing to the directory where the preprocessed data is stored and `config_path` which points to the directory shere the yaml config files are stored which specify which data is considered in the first place. 
We have separate directories for the preprocessed data for every experiment. Thus, to load data for, say lgm and historical experiments, we would call loadData twice with the respective paths as `base_path`-argument.

````julia
import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_to_config_dir = "configs/recipe_configs";

data = sw.loadData(path_to_config_dir, base_path)
````

Then, the function `loadData` has a couple of optional arguments:

- `dir_per_var::Bool` (default: true)
    If true, it is assumed that there is one subdirectory in base_path for each climate variable. If this is not the case, set to false. 

- `is_model_data::Bool` (default: true) for loading observational data, set to false.

- `common_models_across_vars::Bool` (default: false) if true, only the data of those models will be loaded that provide data for all given variables (as retrieved from config files)

- `subset::Dict{String, Vector{String}}` (default: empty)
    
    Specify when you only want to load a subset of the given data. Data is loaded only from files with names that contain any of the mapped values in the dictionary. The following keys are considered (except for the last three, these are fields of the struct `DataID`).

    - `variable`: the short name of the climate variable, e.g. ["tas", "tos"].
    - `statistic`: the statistic/diagnostic that was used when preprocessing the data. The names are generally arbitrary, but need to be identical to those you used in the preprocessing; for instance, we refer to the climatological average as "CLIM", s.t. an example value for this argument is ["CLIM"].
    - `alias`: like for `statistic` the names are arbitrary but must be identical to what you used when preprocessing the data; e.g. ["historical", "historical1"].
    - `exp`: experiment, e.g., ["historical"].
    - `timerange`: timerange; especially for historical data, you may preprocess data for different timeranges,  e.g. ["full", "1980-2014"]. 

    An `alias` refers to a certain `timerange`, e.g. we call the time period from 1951-1980 'historical1'. To load only this data, it thus does not matter whether you set `timerange=["1951-1980"]` or `alias=["historical1"]`.


    - `projects`:  e.g. ["CMIP5", "CMIP6"]. All filenames must contain at least one of the given strings. If not specified and `is_model_data=true`, it is set to ["CMIP"]. If not specified and `is_model_data=false`, it is set to ["ERA5"], which is thus the default observational dataset used.
    - `models`: can refer to models or individual model members, e.g. ["AWI-CM-1-1-MR"]. All filenames must contain at least one of the given strings + "_". The underscore is important since some models have names that are substrings of other models, e.g. "CNRM-CM5" and "CNRM-CM5-C2". 
    - `subdirs`: used to constrain subdirectories to be considered when the data is stored in different subdirectories for each variable (i.e. when `dir_per_var=true`), e.g. ["20241121"] for the case that you only want data from that date and the names of the subdirectories that contain the data include the date.


## How to compute weights

For a full example, see `scripts/run-climwip-simplified.yml`.

Calling sw.computeWeights(model_data, obs_data, config) will compute weights based on the provided data.
The config parameter is of type `ConfigWeights` defined in `src/data-utils.jl`. It holds information 
about the contribution of each combination of statistic/diagnostic and climate variable, once for computing
independence weights and once for computing performance weights. Further parameters from the weighting approach
are specified here (`sigmaD`, `sigmaS`). 

The output of the function `computeWeights` is an object of type `ClimwipWeights` (see `src/data-utils.jl`) which
holds the performance weights for all combinations of statistics/diagnostics and climate variables, for performance as well as
independence weights. 
Further, it contains the overall weights (one for each model, summing up to 1) as well as the performance and independence weights 
for each variable (summed across statistics/diagnostics).
<!-- - `weights_variables:`: For each of 'performance' and 'independence' one value per climate variable considered. These values represent the weight of how much each climate variable influences the generalized distance of a model, which is computed by taking a weighted average across the distances with respect to different variables. Should sum up to 1.  -->

