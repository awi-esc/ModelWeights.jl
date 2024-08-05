include("data.jl")
using LinearAlgebra

@testset "Testset model-model" begin
    data = filter(((k,v),) -> k=="tas", VAR_TO_PREPROC_DATA);
    modelData =  SimilarityWeights.loadPreprocData(data, "CLIM", ["CMIP"]);
    modelDistances = SimilarityWeights.getModelDistances(modelData["tas"])
    
    # make symmetrical matrix
    mat = Array(modelDistances);
    symDistMatrix = mat .+ mat';
    dim = Array(dims(modelData["tas"], :model));
    modelDistancesSym = DimArray(symDistMatrix, (Dim{:model1}(dim), Dim{:model2}(dim)));

    expected = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_tas_CLIM.nc"))["dtas_CLIM"];
    
    # in climwip recipe, the independence model to model distances are not normalized wrt to area weights. 
    # this does not influence the overall weights in the end since the distance matrices are later normalized
    # by the median value across models. So that  the unnormalized values get later normalized by a median
    # that is just a multiple of the median value that the normalized data gets later normalized with. Which then 
    # yields the same result.
    diffFactor = Array(expected)./Array(modelDistancesSym);
    # take care of diagonal values which will be NaN as divided by 0
    d = diffFactor[1,2]
    for i in 1:size(modelDistancesSym)[1]
        diffFactor[i,i] = d
    end
    # the expected and resulting distances should therefore differ only up to a certain factor. 
    @test allequal(round.(diffFactor, digits=NB_DIGITS))
end


@testset "Testset model-obs" begin
    data = filter(((k,v),) -> k=="tas", VAR_TO_PREPROC_DATA);
    models =  SimilarityWeights.loadPreprocData(data, "CLIM", ["CMIP"]);
    observations = SimilarityWeights.loadPreprocData(data, "CLIM", ["ERA5"]);

    modelObsDistances = SimilarityWeights.getModelDataDist(models["tas"], observations["tas"])
    expected = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_tas_CLIM.nc"))["dtas_CLIM"];
    @test round.(Array(modelObsDistances), digits=NB_DIGITS) == round.(Array(expected), digits=NB_DIGITS)
end