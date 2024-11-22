# SimilarityWeights

```@contents
```

This Julia package computes weights for a set of climate models following the approach from Brunner et al (2020). 

## Requirements

#### Preprocessed data
We assume that you have your preprocessed data ready, so that the data from
the different models can be combined, meaning that all models must for instance 
have the same grid.
We do the preprocessing of the data with ESMValTool. A set of recipes to 
download and preprocess some data can be found in our repository [ESMDataPrep](https://github.com/awi-esc/ESMDataPrep).

The structure of the directories where the preprocessed data is stored is
expected by our tool to follow a certain structure. If there is one subdirectory
for each climate variable (this is, for instance, the case when one ESMValTool
recipe is used for a single variable), the structure should adhere to this:

```bash
├── BASE_DIR
│   └── anyname_VAR
│   └── anyname_VAR
│       └── preproc
│       │      └── TASKNAME_VAR
│       │      │       └── VAR_STATISTIC
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      └── TASKNAME_VAR
│       │      │       └── VAR_STATISTIC
│       │      │       └── VAR_STATISTIC
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       └── possibly other output from ESMValTool
....
```

The structure is basically the same if there is not a seperate subdirectory for
each climate variable, except that the BASE_DIR refers to the directory that
immediately contains the preproc-subdirectory: 


```bash
├── BASE_DIR
│       └── preproc
│       │      └── TASKNAME_VAR
│       │      │       └── VAR_STATISTIC
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       │      │       └── VAR_STATISTIC
│       │      └── TASKNAME_VAR
│       │      │       └── VAR_STATISTIC
│       │      │       └── VAR_STATISTIC
│       │      │            └── model1.nc
│       │      │            └── model2.nc
│       │      │            └── ...
│       └── possibly other output from ESMValTool
....
```

#### Config files

Further, to load the data, we need one or more yaml configuration files.
From these, we retrieve which combination of variables, statistics, aliases,
experiments and timeranges will be considered.
   
Since we use ESMValTool to preprocess the data, we can simply use our 
ESMValTool recipes (the fully spelled out version output by ESMValTool, which
is stored in the run-folder and has the name of recipe + "_filled.yml") as 
config files here. Not everything in the recipes is needed for loading the
data. For a minimal example, see `/configs/recipe_configs/recipe_historical_pr_filled.yml` where the unnecessary sections were removed.


## Getting started

For the entire example code, see `scripts/example.jl`.

### How to load data

To load data, you always have to specify the path to the directory that contains the preprocessed data (`base_path`) and the path to the directory that contains the yaml config files (`config_path`).

````julia
import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_to_config_dir = "configs/recipe_configs";

data = sw.loadData(path_to_config_dir, base_path)
````

Then, there are a couple of optional arguments. 

- `dir_per_var::Bool` (default: true)

    If true, it is assumed that there is one subdirectory in base_path for each climate variable. If this is not the case, set to false. 

- `is_model_data::Bool` (default: true) for loading observational data, set to false.

- `common_models_across_vars::Bool` (default: false) if true, only the data of those models will be loaded that provide data for all given variables (as retrieved from config files)

- `subset::Dict{String, Vector{String}}` (default: empty)
    
    Specify when you only want to load a subset of the given data. Data is loaded only from files with names that contain any of the mapped values in the dictionary. The following keys are considered (except for the last three, these are fields of the struct `DataID`).

    - `variable`: refers to short_name, e.g. ["tas", "tos"]
    - `statistic`: e.g. ["CLIM"]. The names of the statistics/diagnostics depend on how you called them when preprocessing the data. 
    - `alias`: e.g. ["historical", "historical1"]
    - `exp`: e.g., ["historical"]
    - `timerange`: e.g. ["full", "1980-2014"]

    An `alias` refers to a certain `timerange`, e.g. we call the time period from 1951-1980 'historical1'. To load only this data, it thus does not matter whether you set `timerange=["1951-1980"]` or `alias=["historical1"]`.

    <br/>

    - `projects`: e.g. ["CMIP5", "CMIP6"]
    - `models`: can refer to models or individual model members, e.g. ["AWI-CM-1-1-MR"]
    - `subdirs`: constrain subdirectories to be considered given that the data is stored in different subdirectories for each variable (argument *dir_per_var* in function *loadData* is true), e.g. ["20241121"]


### How to compute weights

For a full example, see `scripts/run-climwip-simplified.yml`.

Calling ``sw.computeWeights(model_data, obs_data, config)`` will compute weights based on the provided data.
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

#### Defintion of weights: area-weighted RMSE

$d_i^k$ is the area-weighted root mean squared error between the predictions of a model, namely the $k$-th member of model $i$, and the observational data (or for independence weights between pairs of models):

```math
d_i^k = \sqrt{\sum_l w_l (X_l^{a, (i,k)} - X_l^{a, obs})^2}
```

$l$ iterates over (lon, lat)-positions and $w_l$ are the area weights, which depend on the latitudes.
$d_i^k$ is computed for all individual model members, for each diagnostic and variable (denoted as *a*, which may, for instance, refer to the climatological average, CLIM, of near-surface air temperature, tas).



#### Generalized distance

The weights are computed for each model (i.e. *not* on the level of model members). They are based on the **generalized distances**. 
To compute the generalized distances and subsequently the actual weights per model ($w_i$), the different members of each model are first summarized to yield one distance value per model. This is defined as the average distance, $d^\prime_i$ over all distances of $K_i$ members of model $i$ (for a specific combination of diagnostic and variable $a$): 

```math
d_i ^{\prime a} = \frac{\sum_k^{K_i} d_i^k}{K_i}
```

The generalized distances are then defined as follows:

```math
D_i = \sum_a w_a \cdot \frac{d_i^{\prime a}}{\textrm{MEDIAN}(d^a)}
```

Here, the sum iterates over the combination of diagnostics and variables ($a$) and $w_a$ refers to the weights for each combination of diagnostic and variable.
So, $D_i$ is the weighted average over all model distances which are further normalized by the median computed across *all* model members, seperately for each combination of diagnostic and variable. So, $d^a$ in the denominator refers to the distances of all model members for the combination of diagnostic and variable, $a$. 
Note: it doesn't matter if you first average the distances ($d^{\prime a}_i$), to get one value per model and normalize then or normalize first and average then (given that the normalization is in both cases the median across all models and members).

That is, we get a generalized distance value $D_i$ for every model $i$, respectively values $S_{ij}$ for every pair of models $i,j$.

#### Computation of overall weight for model $i$

The generalized distances between models and observations, respectively between model pairs are then combined as follows to yield one weight value for each model $i$:

```math
w_i = \frac{e^{-(\frac{D_i}{\sigma_D})^2}}{1 + \sum_{j \ne i} e^{-\left( \frac{S_{ij}}{\sigma_S} \right)^2}}
```

$\sigma_D$ and $\sigma_S$ are free parameters that Brunner et al. estimated using perfect model tests. For now, we just set them to fix values of 0.5 each.


## References

- Brunner Lukas, Angeline G. Pendergrass, Flavio Lehner, Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming from CMIP6 Projections When Weighting Models by Performance and Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020): 995–1012. https://doi.org/10.5194/esd-11-995-2020.



## Functions

```@autodocs
Modules = [SimilarityWeights]
Order = [:function, :type]
```

## Index

```@index
```

