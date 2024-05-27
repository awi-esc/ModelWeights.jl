using NetCDF

PATH_TO_WORK_DIR = joinpath("recipe_climwip_test_basic_data", "work")

function loadNCdataInDir(path2Dir, climateVar) 
    data = []
    for (root, dirs, files) in walkdir(path2Dir)
        for file in files
            if endswith(file, ".nc") #&& occursin("CMIP", file)
                println(path2Dir * file)
                modelData = NetCDF.ncread(path2Dir * file, climateVar);
                push!(data, modelData);
            end
        end
    end
    # at this point data is a n-element Vector, whose elements each are k-element Vectors
    # data... (splat operator) will unpack the n-element Vector and hcat will concatenate them horitzontally, i.e. 
    # along dimension 2 (i.e. column-wise) so the result is a k x n matrix
    data = hcat(data...);
    return data
end

function _getPath2Data(metric, workDir, varName)
    fn = join([metric, varName, "nc"], "_", ".");
    path2Data = joinpath(workDir, fn);
    print("load data from file: " * path2Data)
    return path2Data
end

function _loadDistMatrix(path2Data, climateVar)
    data = NetCDF.ncread(path2Data, climateVar);
    models = NetCDF.ncread(path2Data, "model_ensemble");
    return data, models
end


function loadIndependentDistMatrix(workDir, varName, climateVar)
    path2Data = _getPath2Data("independence", workDir, varName);
    data, models = _loadDistMatrix(path2Data, climateVar);
    modelRefs = NetCDF.ncread(path2Data, "model_ensemble_reference");
    return data, models, modelRefs
end


function loadPerformanceDistMatrix(workDir, varName, climateVar)
    path2Data = _getPath2Data("performance", workDir, varName);
    data, models = _loadDistMatrix(path2Data, climateVar);
    return data, models
end
