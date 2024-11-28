# SimilarityWeights

```@contents
```

The purpose of this package is twofold: on the one hand, we built it to handle preprocessed Earth System data easily
within Julia and, on the other hand, to compute weights for different sets of Earth System Models following the approach from Brunner et al (2020).

That is, the package provides two types of objects, a data object and a weight object.

## Data Object

The resulting data object has four fields. 
1. `base_path` string that points to the base directory where the preprocessed data is stored.
2. `data_paths` vector of Strings that point to the subdirectories (of base_path) from which the data was loaded.
3. `ids` a list of DataID objects that specify which data was loaded as given by 'variable', 'statistic', experiment, 'exp', 'alias' and 'timerange' where the alias refers to a certain timerange. This shall give a general overview of the loaded data and provide information about the mapping between alias and timerange.

````julia
import SimilarityWeights as sw

# set base_path and config_path

data = sw.loadData(base_path, config_path)

julia> data.ids
2-element Vector{SimilarityWeights.DataID}:
 key=tas_CLIM_lgm variable=tas statistic=CLIM alias=lgm exp=lgm timerange=full 
 key=tos_CLIM_lgm variable=tos statistic=CLIM alias=lgm exp=lgm timerange=full 
````

4. `data` dictionary that contains the actual data. It maps from the keys as given in the DataID objects (a string concatenating variable, statistic and alias with underscores) to the respective data which each is a DimensionalData.DimArray.


## Weight Object

````julia
import SimilarityWeights as sw

# load model and observational data

config_weights = sw.ConfigWeights(
    performance = Dict("tas_CLIM" => 1, "tos_CLIM" => 1),
    independence = Dict("tas_CLIM" => 1, "tos_CLIM" => 1),
    sigma_independence = 0.5,
    sigma_performance = 0.5,
    ref_period = "historical", # can be alias or timerange
    target_dir = "/weights/"
);

weights = sw.computeWeights(model_data_historical, obs_data, config_weights);
````

The `ClimwipWeights` object has several fields that each store a
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
Modules = [SimilarityWeights]
Order = [:function, :type]
```

## Index

```@index
```

