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
