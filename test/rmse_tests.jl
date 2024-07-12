include("data.jl")
using NCDatasets
using LinearAlgebra

NB_DIGITS = 5

@testset "Testset model-model" begin
    modelData =  SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas"], "CLIM", ["CMIP"]);
    modelDistances = SimilarityWeights.getModelDistances(modelData["tas"])
    
    expected = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_tas_CLIM.nc"))["dtas_CLIM"];
    
    # in climwip recipe, the independence model to model distances are not normalized wrt to area weights. 
    # this does not influence the overall weights in the end since the distance matrices are later normalized
    # by the median value across models. So that  the unnormalized values get later normalized by a median
    # that is just a multiple of the median value that the normalized data gets later normalized with. Which then 
    # yields the same result.
    diffFactor = Array(expected)./Array(modelDistances);
    # take care of diagonal values which will be NaN as divided by 0
    d = diffFactor[1,2]
    for i in 1:size(modelDistances)[1]
        diffFactor[i,i] = d
    end
    # the expected and resulting distances should therefore differ only up to a certain factor. 
    @test allequal(round.(diffFactor, digits=NB_DIGITS))
end


@testset "Testset model-obs" begin
    models =  SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas"], "CLIM", ["CMIP"]);
    observations = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas"], "CLIM", ["ERA5"]);

    modelObsDistances = SimilarityWeights.getModelDataDist(models["tas"], observations["tas"])
    expected = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_tas_CLIM.nc"))["dtas_CLIM"];
    @test round.(Array(modelObsDistances), digits=NB_DIGITS) == round.(Array(expected), digits=NB_DIGITS)
end