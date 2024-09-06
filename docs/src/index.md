# SimilarityWeights

```@contents
```

This Julia package computes weights for a set of climate models following the approach
from Brunner et al (2020). 

## Getting started

Here's an example (in Julia) of how to get weights for a set of models defined in the configuration 
file stored in configs/: 

````julia
using SimilarityWeights

path_config = "configs/example_historical_albedo.yml"
config = SimilarityWeights.validateConfig(path_config);

weights, avgs = SimilarityWeights.runWeights(config);
````

#### Configuration

The assumed directory structure of the data is as follows (upper case words are variables 
that need to be replaced): 

```bash

├── EXPERIMENT
│   └── [PREFIX_VAR_FOLDERS_]VAR
│       ├── preproc
│       │   └── climatology_full
│       │   │   └── VAR
│       │   ├── climatology_historical1
│       │   │   └── VAR
│       │   ├── climatology_historical2
│       │   │   └── VAR
│       │   └── climatology_historical3
│       │       └── VAR
│       │     
│       └── run
│
└── EXPERIMENT
│   └── ...
....

```
Note that the 'run'-directory will be there if the data was loaded with our ESMValTool recipes, but it doesn't contain any data that we'll need. So no need for this directory if the data was loaded differently. 


- `base_path:`  Path to the directory that contains the preprocessed data from ESMValTool (not necessarily from ESMValTool, but the underlying structure must be the same)

- `experiment:` Name of experiment for which models were run, e.g. 'historical' or 'midHolocene'.

- `prefix_var_folders:` (optional) The prefix of the directories for each climate variable. 

- `target_dir:` Path to the directory where the computed data will be stored.

- `variables:` List of climate variables which will be considered.

- `name_ref_period:` If the data is loaded with our ESMValTool recipes, the name of the reference period is, for now, one of 'historical1', 'historical2', 'historical3' (see below).

- `name_full_period:` If the data is loaded with our ESMValTool recipes, the name with which the entire referenced period of the respective experiment is referred to, is set to 'full'.

- `models_project_name:` Either 'CMIP6' or 'CMIP5'. We focus on 'CMIP6'.

- `obs_data_name:` If the data is loaded with our ESMValTool recipes, for now this is set to 'ERA5'. 

- `weight_contributions:` For now: one value for performance, one for independence. Should sum up to 1. This is how much each of the two weight-types is taken into account.

- `weights_variables:`: For each of 'performance' and 'independence' one value per climate variable considered. These values represent the weight of how much each climate variable influences the generalized distance of a model, which is computed by taking a weighted average across the distances with respect to different variables. Should sum up to 1. 

## Computation of weights

#### Distances: area-weighted RMSE
``d_i``is the area-weighted root mean squared error between model predictions and observational data.
This values is computed for every model (including different ensemble members)[^1], diagnostic and variable (summarized as DIAG):

```math
d_i = \sqrt{\sum_l w_l (X_l^{DIAG,MODEL:i} - X_l^{DIAG, Obs})^2}
```
[^1]: Note that I write **model** in lower-case when I refer to all models, including different ensemble members. When I refer to all models where each is a summary of its ensemble members, I write **Model** in upper-case.

#### Generalized distance

```math
D_i = \sum_a \frac{w_a \cdot d_i}{\textrm{MEDIAN}(d_i^a)}
```

Here, the sum iterates over the combination of diagnostics and variables and represents the **generalized distance** for a model ``i``.
These distances are normalized by their median values.
For example for the diagnostic climatological average and variable sea surface temperature, we compute for every model 
the distance ``d_i`` and normalize it by the median distance across all models. 
Then to compute the generalized distance of a model ``i``, we take the weighted average over each combination of diagnostic and variable using weights ``w_a``.  
That is, we get a generalized distance value ``D_i`` for every model, including different ensemble members.

#### Overall weight per Model

```math
w_i = \frac{e^{-(\frac{D_i}{\sigma_D})^2}}{1 + \sum_{j \ne i}^{M} e^{-\left( \frac{S_{ij}}{\sigma_S} \right)^2}}
```

To compute the actual weight per Model, ``w_i``, the different ensemble members are first summarized to yield one distance value per Model, which is the average distance, ``d^\prime_i`` over all distances of ``K_i``ensemble members of Model ``i``: 

```math
d_i ^\prime = \frac{\sum_k^{K_i} d_i^k}{K_i}
```


#### Reference data

We use three different periods within the historical period as reference to compute the performance weights: 

- historical1: 1951 - 1980
- historical2: 1961 - 1990
- historical3: 1990 - 2014


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

