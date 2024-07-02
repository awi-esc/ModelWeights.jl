using NetCDF

PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "work")
PATH_TO_PREPROC_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "preproc")
"""
    loadNCdataInDir(path2Dir, climateVar[, dataIncluded=[], addHorizontal=true])

Loads all .nc files inside the directory loacated at 'path2Dir' into a single Vector.

If length(dataIncluded) != 0 only those .nc files are considered whose filenames contain all elements of dataIncluded, 
e.g. if dataIncluded=['ERA5'] only ERA5 data will be included (files with ERA5 in their filename).
"""
function loadNCdataInDir(path2Dir, climateVar, dataIncluded = []) 
    data = []
    sources = []
    for (root, dirs, files) in walkdir(path2Dir)
        for file in files
            if endswith(file, ".nc")
                if length(dataIncluded) != 0
                    if all([occursin(name, file) for name in dataIncluded])
                        print(file)
                        push!(sources, file)
                        modelData = NetCDF.ncread(joinpath(path2Dir, file), climateVar);
                        push!(data, modelData);
                        println("added model data: " * string(size(modelData)))
                    end
                else
                    print(file)
                    push!(sources, file)
                    modelData = NetCDF.ncread(joinpath(path2Dir, file), climateVar);
                    push!(data, modelData);
                    println("added model data: " * string(size(modelData)))
                end
            end
        end
    end
    return data, sources
end

function _getPath2Data(metric, workDir, varName)
    fn = join([metric, varName, "nc"], "_", ".");
    path2Data = joinpath(workDir, fn);
    println("load data from file: " * path2Data)
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


function loadPerformanceMetric(workDir, varName, climateVar)
    path2Data = _getPath2Data("performance", workDir, varName);
    data, models = _loadDistMatrix(path2Data, climateVar);
    return data, models
end
