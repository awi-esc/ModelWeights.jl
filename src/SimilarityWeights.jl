module SimilarityWeights

include("data-utils.jl")
include("plot-utils.jl")
include("weights.jl")
include("plot-weights.jl")
include("plot-data.jl")
include("main-functions.jl")

export loadPreprocData, getIndependenceWeights, getPerformanceWeights, getOverallWeights, runWeights

end
