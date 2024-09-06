# SimilarityWeights

```@contents
```

This Julia package computes weights for a set of climate models following the approach
from Brunner et al (2020). 

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

