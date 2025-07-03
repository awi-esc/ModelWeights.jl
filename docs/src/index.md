# Usage

```@contents
```

The purpose of this package is twofold: on the one hand, we built it to handle preprocessed Earth System data easily
within Julia and, on the other hand, to compute weights for different sets of Earth System Models following the approach from Brunner et al (2020) and 
adding the possibility to use and combine different performance metrics.

## Data Representation

The basic representation of data is in form of a dictionary that simply maps from ids 
(Strings) to the data itself stored as YAXArrays.
We name this type alias a `DataMap`. 

The metadata of each dataset (``datamap[id].properties``) contains the 
metadata loaded from the original files as well as additional information if it is provided
in the optional argument `meta_data`.


## Weights

<!-- ````julia
import ModelWeights as mw
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
- `w` is the overall weight vector of length n. -->


## Application of weights

Add example of how weights are applied to data


## Functions

```@autodocs
Modules = [ModelWeights.Data, ModelWeights.Plots, ModelWeights.Timeseries, ModelWeights.Weights]
Order   = [:function, :type]
```

## Index

```@index
```

