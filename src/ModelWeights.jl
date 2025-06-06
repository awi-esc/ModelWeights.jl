module ModelWeights

# include submodules
include("Data.jl")
include("Timeseries.jl")
include("Weights.jl")
include("Plots.jl")

# make submodules available in the scope of main module (ModelWeights)
using .Data  
using .Timeseries
using .Weights
using .Plots

export loadDataFromYAML, loadDataFromESMValToolRecipes, climwipWeights


export joinDataMaps
export writeDataToDisk, readDataFromDisk

export getUncertaintyRanges,
    globalMeans,
    addAnomalies!,
    approxAreaWeights,
    anomaliesGM!,
    addAnomaliesGM!
export addLinearTrend!, linearTrend
export setToSummarizedMembers!, summarizeEnsembleMembersVector
export getLandMask, getOceanMask, addMasks!, subsetModelData, alignPhysics


export weightedAvg, applyWeights, getModelLogLikelihoods
export equalWeights, distributeWeightsAcrossMembers
export writeWeightsToDisk, saveWeightsAsNCFile



end # module ModelWeights
