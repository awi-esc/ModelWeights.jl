include("data.jl")

CLIMATE_VARS = ["tas"];
DIAGNOSTIC = "CLIM";
MODEL_DATA =  SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["CMIP"]);

@testset "Testset performance weights" begin    
    obsData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["ERA5"]);
    weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    weights = SimilarityWeights.getPerformanceWeights(MODEL_DATA, obsData, weightsVars);

    
    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
    expected =  Array(ds["overall_mean"]);
    
    @test round.(Array(expected), digits=NB_DIGITS) == round.(Array(weights), digits=NB_DIGITS);
end


@testset "Testset independence weights" begin
    weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    weights = SimilarityWeights.getIndependenceWeights(MODEL_DATA, weightsVars);
    weights = round.(Array(weights), digits=NB_DIGITS);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
    expected = round.(ds["overall_mean"][:,:], digits=NB_DIGITS);

    @test !any((x) -> x!=1, expected .== weights)
end