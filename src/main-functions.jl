import YAML
using DimensionalData


"""
    computeWeights(model_data::Data, obs_data::Data, config::ConfigWeights)

Compute weight for each model in multi-model ensemble according to approach
from Brunner, Lukas, Angeline G. Pendergrass, Flavio Lehner,
Anna L. Merrifield, Ruth Lorenz, and Reto Knutti. “Reduced Global Warming
from CMIP6 Projections When Weighting Models by Performance and
Independence.” Earth System Dynamics 11, no. 4 (November 13, 2020):
995–1012. https://doi.org/10.5194/esd-11-995-2020.

# Arguments:
- `model_data`: Models for data for computing independence and performance weights.
- `obs_data`: Observational data for computing performance weights.
- `config`: parameters specifiying the relative contributions of each 
combination of variable and diagnostic.
"""
function computeWeights(
    model_data::Data, obs_data::Data, config::ConfigWeights
)
    weights_perform = normalizeWeightsVariables(config.performance)  
    weights_indep = normalizeWeightsVariables(config.independence)
    
    # sanity checks for input arguments
    keys_weights_perform = allcombinations(dims(weights_perform, :variable), dims(weights_perform, :diagnostic))
    keys_weights_indep = allcombinations(dims(weights_indep, :variable), dims(weights_indep, :diagnostic))


    ref_period_alias = getRefPeriodAsAlias(config.ref_period)
    msg =  x -> "For computation of $x weights: Make sure that data is provided 
    for the given reference period ($(config.ref_period)) and combination of 
        variables and diagnostic for which (balance) weights were specified!"
    if !isValidDataAndWeightInput(model_data, keys_weights_perform, ref_period_alias)
        throw(ArgumentError(msg("performance")))
    end
    if !isValidDataAndWeightInput(obs_data, keys_weights_perform, ref_period_alias)
        throw(ArgumentError(msg("performance")))
    end
    if !isValidDataAndWeightInput(model_data, keys_weights_indep, ref_period_alias)
        throw(ArgumentError(msg("independence")))
    end

    dists_perform_all = computeDistancesAllDiagnostics(
        model_data, obs_data, config.performance, ref_period_alias, true
    )
    dists_indep_all = computeDistancesAllDiagnostics(
        model_data, obs_data, config.independence, ref_period_alias, false
    )
    Di = computeGeneralizedDistances(dists_perform_all, weights_perform, true)
    Sij = computeGeneralizedDistances(dists_indep_all, weights_indep, false)

    performances = performanceParts(Di, config.sigma_performance)
    independences = independenceParts(Sij, config.sigma_independence)
    weights = performances ./ independences;
    weights = weights ./ sum(weights);
    setRefPeriodInWeightsMetadata!(weights.metadata, config.ref_period, ref_period_alias)
    
    wI = independences ./ sum(independences)
    wP = performances ./ sum(performances)
    setRefPeriodInWeightsMetadata!(wP.metadata, config.ref_period, ref_period_alias)


    climwip_weights =  ClimwipWeights(
        performance_distances = dists_perform_all,
        independence_distances = dists_indep_all, 
        Di = Di,
        Sij = Sij,
        wP = wP,
        wI = wI,
        w =  weights
        #overall = (wP./wI)./sum(wP./wI), # just for sanity check
    )
    logWeights(weights.metadata);
    if !isempty(config.target_dir)
        saveWeights(climwip_weights, config.target_dir)
    end
    return climwip_weights
end
