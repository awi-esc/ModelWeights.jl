module SimilarityWeights

include("main-functions.jl")
include("weights.jl")
include("plot-utils.jl")
include("plot-weights.jl")
include("plot-data.jl")
include("data-utils.jl")

export loadPreprocData, getIndependenceWeights, getPerformanceWeights, getOverallWeights, runWeights

end
