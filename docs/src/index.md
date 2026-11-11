# ModelWeights.jl

`ModelWeights.jl` was built to facilitate working with climate model ensembles such as CMIP, with a particular focus on computing model weights.

Among other things, it provides functions to compute model weights based on independence and performance metrics, following the approach of Brunner et al. (2020).

To use the provided functionalities, it is assumed that the data is already available in a preprocessed state, in particular, it must be defined on the same lon-lat grid.


## Package features

- load data for different variables/experiments into convenient data format
- filter data to get only models with specific variables for various experiments
- compute independence weights based on model similarity
- compute performance weights based on historical period or based on custom function
- combine independence- and performance weights

## Content

```@contents
```






