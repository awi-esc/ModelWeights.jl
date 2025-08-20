# Model weights

The weighting method that we use was developed by Knutti et al. (2017) that is based on work from Sanderson et al. (2015) (see [References](@ref)).

## Basic computation of weights
 
The overall weight for model $i$ is computed based on its performance- and independenct weights, $w^P_i$ and $w^I_i$:

```math
w_i = a_0 \cdot w^P_i \cdot w^I_i 
```

The performance and independence weights are defined as follows: 

```math
w^{P}_i = a_1 \cdot e^{-(\frac{D_i}{\sigma_D})^2}
```

```math
w^{I}_{i} = a_2 \cdot \frac{1}{\sum_j e^{-\left( \frac{S_{i,j}}{\sigma_S} \right)^2}} 
```

$S_{i,j}$ refers to the generalized distances between models $i$ and $j$ and $D_i$ refers to the generalized distances between predictions of model $i$ and observatioanl/reanalysis data
$\sigma_S$ and $\sigma_D$ are hyperparameters that specify the strength of the independence- and performance weighting respectively and $a_0, a_1, a_2$ are normalization constants.


<!--- 
## Combination of performance weights

One may want to measure performance based on data other than data from the historical period and potentially combine different performance measurements.
We provide the possibility to compute performance weights based on a custom function and combine them into one overall performance weight:

```math
w_i^P = a_3 \prod_k \cdot w^{P_k}_i 
```




## Area-weighted RMSE

```math
d_i^k = \sqrt{\sum_l w_l (X_l^{a, (i,k)} - X_l^{a, obs})^2}
```

The variable $d_i^k$ is the area-weighted root mean squared error between the predictions of a model, namely the $k$-th member of model $i$, and the observational data (or for independence weights between pairs of models):

Variable $l$ iterates over (lon, lat)-positions and $w_l$ are the area weights, which account for the different grid cell areas depending on the latitudes.
$d_i^k$ is computed for all individual model members, for each diagnostic and variable (denoted as *a*, which may, for instance, refer to the climatological average, CLIM, of near-surface air temperature, tas).


## Generalized distances

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

## Computation of overall weight for model $i$

The generalized distances between models and observations, respectively between model pairs are then combined as follows to yield one weight value for each model $i$ ($w_i$'s are normalized by $\sum_i w_i$):

```math
w_i = \frac{e^{-(\frac{D_i}{\sigma_D})^2}}{1 + \sum_{j \ne i} e^{-\left( \frac{S_{ij}}{\sigma_S} \right)^2}}
```

The parameters, $\sigma_D$ and $\sigma_S$, are free parameters that Brunner et al. estimated using perfect model tests.

## Computation of performance/independence weights only

When computing only performance weights, we assume that all models are equally dependent:

```math
w^{P}_i = \frac{e^{-(\frac{D_i}{\sigma_D})^2}}{N}
```

When computing only independence weights, we assume that all models perform equally well.

```math
w^{I}_i = \frac{1}{1 + \sum_{j \ne i} e^{-\left( \frac{S_{ij}}{\sigma_S} \right)^2}}
```

---> 
