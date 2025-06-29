# Usage

```@contents
```

The purpose of this package is twofold: on the one hand, we built it to handle preprocessed Earth System data easily
within Julia and, on the other hand, to compute weights for different sets of Earth System Models following the approach from Brunner et al (2020).

The main objects that the package provides are: `DataMap` and two weight 
objects, `ConfigWeights` and `Weights`.


## DataMap Object

Each Data object refers to some summary data (e.g. climatological average) 
of a particular experiment for the full time period or part of it.

`DataMap` is a type alias for a dictionary mapping from Strings to YAXArrays.
The metadata of each dataset (``datamap[id].properties``) contains the 
metadata loaded from the original files as well as some additional information 
that we add: 
- `_variable`
- `_statistic`
- `_alias`
- `_timerange`
- `_paths`
- `_id`: concatenates `variable`, `statistic` and `alias` into a single String 
separated by underscores


## Weight Objects

### ConfigWeights

````julia
import ModelWeights as mw

# load model and observational data

config_weights = mw.ConfigWeights(
    performance = Dict("tas_CLIM" => 1, "tos_CLIM" => 1),
    independence = Dict("tas_CLIM" => 1, "tos_CLIM" => 1),
    sigma_performance = 0.5,
    sigma_independence = 0.5,
    alias_ref_perform_weights = "historical",
    alias_ref_indep_weights = "historical",
    target_path = ""
);
````

### Weights

````julia
weights = mw.computeWeights(model_data_historical, obs_data, config_weights);
````

The `Weights` object has several fields that each store a
DimensionalData.DimArray. Some of them refer to normalized weights, others refer 
to the distances on the basis of which the weights were computed.

- `performance_distances` and `independence_distances` contain the distances
between models and data (performance) or between models and models (independence) 
for every combination of variable and diagnostic/statistic.
- `Di` and `Sij` are the generalized distances, i.e. weighted average of all
distances across variables and diagnostic. So, $D_i$ is a vector with size 
1xn where n is the number of models (not on level of members!) and $S_{ij}$ is 
a matrix of size n x n. 

- `wP` and `wI` respectively store the normalized performance/independence weights, 
both have size 1xn (n: number of models).
- `w` is the overall weight vector of length n.


## Application of weights

TODO: add example of how weights are applied to data


## Functions

```@autodocs
Modules = [ModelWeights.Data, ModelWeights.Plots, ModelWeights.Timeseries, ModelWeights.Weights]
```

## Index

```@index
```

