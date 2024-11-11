module SimilarityWeights

const MODEL_MEMBER_DELIM = "#"

include("data-utils.jl")
include("data-functions.jl")
include("plot-utils.jl")
include("weights.jl")
include("plot-weights.jl")
include("plot-data.jl")
include("main-functions.jl")

export loadData, loadPreprocData, computeWeights, getCommonModelsAcrossVars

end
