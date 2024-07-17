include("data.jl")

# CLIMATE_VARS = ["tas"];
# DIAGNOSTIC = "CLIM";
# MODEL_DATA =  SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["CMIP"]);
# OBS_DATA = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, CLIMATE_VARS, DIAGNOSTIC, ["ERA5"]);

# VARIABLE_CONTRIBUTIONS = Dict(
#     "independence" => Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0), 
#     "performance" => Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1)
# );

@testset "Testset performance weights" begin    
    # weights = SimilarityWeights.getPerformanceWeights(MODEL_DATA, OBS_DATA, VARIABLE_CONTRIBUTIONS["performance"]);
    weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);
    obsData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["ERA5"]);
    
    weightsByVars = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
    weights = SimilarityWeights.summarizeWeightsAcrossVars(weightsByVars);
    weightsRounded = round.(Array(weights), digits=NB_DIGITS);

    #weightsEnsemble = SimilarityWeights.averageEnsembleVector(weights);
    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "performance_overall_mean.nc"));
    expected =  Array(ds["overall_mean"]);
    expectedRounded = round.(expected, digits=NB_DIGITS);

    @test all((x)-> x==1, expectedRounded .== weightsRounded);
end


@testset "Testset independence weights" begin
    # weights = SimilarityWeights.getIndependenceWeights(MODEL_DATA, VARIABLE_CONTRIBUTIONS["independence"]);
    weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);

    weightsByVars = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
    weights = SimilarityWeights.summarizeWeightsAcrossVars(weightsByVars);
    weights = round.(Array(weights), digits=NB_DIGITS);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "independence_overall_mean.nc"), "r");
    expected = round.(ds["overall_mean"][:,:], digits=NB_DIGITS);

    @test all((x)-> x==1, expected .== weights);
end


@testset "Testset combined weights" begin
    nb_digits = 2;
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);
    weightsVars = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    wI_vars = SimilarityWeights.getIndependenceWeights(modelData, weightsVars);
    wI = SimilarityWeights.summarizeWeightsAcrossVars(wI_vars);


    obsData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["ERA5"]);
    weightsVars = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    wP_vars = SimilarityWeights.getPerformanceWeights(modelData, obsData, weightsVars);
    wP = SimilarityWeights.summarizeWeightsAcrossVars(wP_vars);
    weights = SimilarityWeights.combineWeights(wP, wI, 0.5, 0.5)

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"));
    expected = Array(ds["weight"]);

    #last four entries are from one ensemble 
    weightCCSM4 = weights[4]/4;
    result = [Array(weights[1:3]); repeat([weightCCSM4], 4)]

    models = [Array(dims(weights, :model)); repeat(["CCSM4"], 3)];
    weightsEnsemble = DimArray(result, (Dim{:model}(models)))

    @test all((x)-> x==1, round.(expected[1:3], digits=nb_digits) .== round.(weightsEnsemble[1:3], digits=nb_digits));
end


@testset "Testset function getWeights" begin
    nb_digits = 2;
    modelData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["CMIP"]);
    obsData = SimilarityWeights.loadPreprocData(PATH_TO_PREPROC_DIR, ["tas", "pr", "psl"], "CLIM", ["ERA5"]);
    
    weightsVarsPerform = Dict{String, Number}("tas" => 1, "pr" => 2, "psl" => 1); 
    weightsVarsIndep = Dict{String, Number}("tas" => 0.5, "pr" => 0.25, "psl" => 0); 
    weights = SimilarityWeights.getWeights(modelData, obsData, 0.5, 0.5, weightsVarsPerform, weightsVarsIndep);

    ds = NCDataset(joinpath(PATH_TO_WORK_DIR, "calculate_weights_climwip", "climwip", "weights.nc"));
    expected = Array(ds["weight"]);

    #last four entries are from one ensemble 
    weightCCSM4 = weights[4]/4;
    result = [Array(weights[1:3]); repeat([weightCCSM4], 4)]

    models = [Array(dims(weights, :model)); repeat(["CCSM4"], 3)];
    weightsEnsemble = DimArray(result, (Dim{:model}(models)))

    @test all((x)-> x==1, round.(expected[1:3], digits=nb_digits) .== round.(weightsEnsemble[1:3], digits=nb_digits));
end


