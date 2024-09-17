using NCDatasets

PATH_TO_PREPROC_DIR = joinpath(@__DIR__, "..", "reproduce-climwip-figs", 
    "recipe_climwip_test_basic_data", "preproc", "calculate_weights_climwip");

PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "reproduce-climwip-figs", 
    "recipe_climwip_test_basic_data", "work");

# make dictionary from variabes to be tested to paths 
VAR_TO_PREPROC_DATA = Dict{String, Dict{String, String}}();
VAR_TO_PREPROC_DATA["CLIM"] = Dict{String, String}();
for var in ["tas", "pr", "psl"]
    VAR_TO_PREPROC_DATA["CLIM"][var] = joinpath(PATH_TO_PREPROC_DIR, var * "_CLIM"); 
end

# Round results to this nb of digits to compare with actual data
NB_DIGITS = 6;