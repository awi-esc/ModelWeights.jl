module ModelWeights

# using Logging
# DEBUG_LOGGER = ConsoleLogger(stderr, Logging.Debug)



const MODEL_MEMBER_DELIM = "#"
const MODEL_NAME_FIXES = Dict(
    "FGOALS_g2" => "FGOALS-g2",
    "ACCESS1.3" => "ACCESS1-3"
)

include("helper-functions.jl")
include("data-utils.jl")
include("diagnostics.jl")
include("data-functions.jl")
include("plot-utils.jl")
include("weights.jl")
include("plot-weights.jl")
include("plot-data.jl")
include("main-functions.jl")



export loadDataFromYAML, loadDataFromESMValToolConfigs, computeWeights


export joinDataMaps
export writeDataToDisk, readDataFromDisk

export getUncertaintyRanges, computeGlobalMeans, addAnomalies!, approxAreaWeights, computeAnomaliesGM!, addAnomaliesGM!
export addLinearTrend!, computeLinearTrend
export averageEnsembleMembers!, summarizeEnsembleMembersVector
export getLandMask, getOceanMask, addMasks!, subsetModelData, alignPhysics


export computeWeightedAvg, applyWeights, getModelLogLikelihoods
export makeEqualWeights, distributeWeightsAcrossMembers
export writeWeightsToDisk, saveWeightsAsNCFile


end # module ModelWeights
