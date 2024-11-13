# SimilarityWeights

```@contents
```

This Julia package computes weights for a set of climate models following the approach from Brunner et al (2020). 

## Prerequisite Data structure

We assume that you have your preprocessed data ready, so that the data from
the different models can be combined, meaning that all models must for instance 
have the same grid.
We do the preprocessing of the data with ESMValTool. A set of recipes to 
download and preprocess some data can be found in our repository [ESMDataPrep](https://github.com/awi-esc/ESMDataPrep).

The structure of the directroies where the preprocessed data is stored is
expected by our tool to follow a certain structure. If there is one subdirectory
for each climate variable, the structure should adhere to this:

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


Further, to load the data, we need one or more yaml configuration files.
We use ESMValTool to preprocess the data and simply use our ESMValTool recipes
as config files here. Not everything in the recipes is needed for loading the
data. For a minimal example, see `/configs/recipe_configs/recipe_historical_pr_filled.yml` 
where the unnecessary sections were removed.



## Getting started

For the entire example code, see `scripts/example.jl`.

#### How to load data

To load data, you always have to specify the path to the directory that contains the
yaml config files `path_to_config_dir` and the path to the directory that 
contains the prerpocessed data `base_path` (see above for details on structure).

````julia
import SimilarityWeights as sw

base_path = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical";
path_to_config_dir = "configs/recipe_configs";

data = sw.loadData(path_to_config_dir, base_path)
````

Here it is assumed that there is one subdirectory in base_path for each climate variable.
If this is not the case (e.g. see example in `scripts/run-climwip-simplified.yml`), 
set the parameter `dir_per_var=false` when calling `sw.loadData`. 

You might not want to load all data, but only a subset of it. In that case, 
you can specify a `DataConstraint` (struct defined in src/data-utils.jl) like so:

````julia
dc = sw.DataConstraint(
    variables=["tas", "pr"], 
    tasks=["historical1"],
    commonModelsAcrossVars=true
);
````

You can specify which variables to load and specify whether only those models should be
loaded that provide data for all variables (parameter `commonModelsAcrossVars`), which 
is by default set to false.
Further, you can filter the timeranges (`timeranges`, e.g. set to `"1951-1980"`) and 
which tasks (`tasks`) to load. A task usually refers to a certain time period, e.g. 
we call the time period from 1951-1980 'historical1'. 
So in this case it does not matter whether you provied `timerange=1951-1980` or 
`task=historical1`.
Then you can filter the statistics to load (e.g., `statistics=["CLIM"]` or `statistics=["ANOM"]`). 
The names of the statistics/diagnostics depend on how you called them when preprocessing the data. 


#### How to compute weights

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

#### Defintion of the weights: area-weighted RMSE
``d_i``is the area-weighted root mean squared error between model predictions and observational data.
This values is computed for every model (including different ensemble members)[^1], diagnostic and variable (summarized as DIAG):

```math
d_i = \sqrt{\sum_l w_l (X_l^{DIAG,model_i} - X_l^{DIAG, Obs})^2}
```
[^1]: Note that I write **model** in lower-case when I refer to all models, including different ensemble members. When I refer to all models where each is a summary of its ensemble members, I write **Model** in upper-case or refer to it as **ensemble**.

#### Generalized distance

The weights are computed for each Model/ensemble. They are based on the **generalized distances**. 
To compute the generalized distances and subsequently the actual weights per Model/ensemble, ``w_i``, the different ensemble members are first summarized to yield one distance value per Model. This is defined as the average distance, ``d^\prime_i`` over all distances of ``K_i``ensemble members of Model ``i``: 

```math
d_i ^\prime = \frac{\sum_k^{K_i} d_i^k}{K_i}
```

The generalized distances are then defined as follows:

```math
D_i = \sum_a \frac{w_a \cdot d_i^{\prime a}}{\textrm{MEDIAN}(d_i^a)}
```

Here, the sum iterates over the combination of diagnostics and variables (*a*) and `w_a` refers to the weights for each combination of diagnostic and variable.
So, it computes the weighted average over all Model/ensemble distances which are further normalized by the median, computed across *all* models (on level of ensemble members), for each combination of diagnostic and variable that is used to compute D_i or, respectively, S_{ij}.
Note: it doesn't matter if you first average the distances, to get one value per Model/ensemble and normalize then or normalize fist and average then (given that the normalization is in both cases the median across all models/ensemble members).

That is, we get a generalized distance value ``D_i``, or respectively ``S_ij`` for every Model/ensemble (or pair i,j of such).

#### Overall weight per Model

```math
w_i = \frac{e^{-(\frac{D_i}{\sigma_D})^2}}{1 + \sum_{j \ne i}^{M} e^{-\left( \frac{S_{ij}}{\sigma_S} \right)^2}}
```



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

