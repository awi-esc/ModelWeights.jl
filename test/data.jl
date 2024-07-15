using NCDatasets

PATH_TO_PREPROC_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "preproc", "calculate_weights_climwip");
PATH_TO_WORK_DIR = joinpath(@__DIR__, "..", "recipe_climwip_test_basic_data", "work");

# Round results to this nb of digits to compare with actual data
NB_DIGITS = 6;