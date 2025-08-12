# ModelWeights.jl

A package to compute weights (based on independence- and performance) for ensembles of Earth System Models (ESMs), especially CMIP models, including convenience functions for handling climate model data.

## Package features

- only requirement: the data must already be preprocessed, i.e. defined on the same grid
- load data for different variables/experiments into convenient data format
- filter data to get only models with specific variables for various experiments
- compute independence weights based on model similarity
- compute performance weights based on historical period or based on custom function
- combine independence- and performance weights


```@contents
```






