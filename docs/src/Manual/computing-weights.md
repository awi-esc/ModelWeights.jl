# Computing weights
For a full example, see `scripts/run-climwip-simplified.yml`.

Calling `mw.computeWeights(dists_indep, dists_perform, config)` will compute weights based on the distances for the independence weights and the distances for the performance weights.
The config parameter is of type `ConfigWeights` defined in `src/data-utils.jl`. It holds information about the contribution of each combination of statistic/diagnostic and climate variable, once for computing independence weights and once for computing performance weights. Further parameters concerning the weighting are specified in the ConfigWeights struct, such as the hyperparameters,  `sigma_performance` and `sigma_independence`. 


```
The output of the function `computeWeights` is an object of type `Weights` (see `src/data-utils.jl`) which contains the independence and performance weights seperately as well as the overall weights (`wI`, `wP` and `w`), each of which respectively sum up to 1. Further it contains the following data:
- For all combinations of statistics/diagnostics and climate variables, the distances used for computing performance as well as independence weights (`performance_distances`, `independence_distances`). 
- For every model, the generalized distances for performance (`Di`) and independence (`Sij`), (distances summed across statistics/diagnostics + variables).
- Model weights on the basis of members, i.e. the weight for every model is distributed evenly across its members, so this gives one weight per member.
- the configuration in form of the `ConfigWeights` object used to compute the weights.


```@raw html
<!-- `weights_variables:`: For each of 'performance' and 'independence' one value per climate variable considered. These values represent the weight of how much each climate variable influences the generalized distance of a model, which is computed by taking a weighted average across the distances with respect to different variables. Should sum up to 1.  -->
```



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

s