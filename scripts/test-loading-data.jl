import ModelWeights as mw
import ModelWeights.Data as mwd

using BenchmarkTools
using DimensionalData
using Profile, ProfileSVG
using YAXArrays

data_dir = "./data/"
base_dir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool"
paths_data = joinpath.(
    base_dir,
    [
     "historical/recipe_cmip6_historical_tas_timeseries_20250228_081213/preproc/historical/tas_CLIM-ann",
     #"ssp585/recipe_cmip6_ssp585_tas_timeseries_20250226_074213/preproc/ssp585/tas_CLIM-ann",
     #"historical/recipe_cmip6_historical_psl_timeseries_20250313_130922/preproc/historical/psl_CLIM-ann"
    ]
)
data_ids = [
    "tas_annual_historical", 
    #"tas_annual_ssp585", 
    #"psl_annual_historical"
]; 


# ---------------------------- Preview Model Data ---------------------------------------- #
# CMIP6 data
preview = mwd.previewDataMap(
    paths_data,
    data_ids,
    level = :model,
    filename_format = :esmvaltool_cmip6,
    constraint = Dict(:models => ["AWI-CM-1-1-MR"]),
    #constraint_ts = (start_year = 1950, end_year = 2014) # for this filename format constraint_ts no influence on preview data
)

# CMIP data 
preview = mwd.previewDataMap(
    ["/albedo/work/projects/p_pool_clim_data/CMIP6/CMIP/AWI/AWI-CM-1-1-MR/historical/r1i1p1f1/Amon/tas/gn/v20200720"], 
    ["tas_annual_historical"],
    level = :model,
    dtype = "cmip",
    filename_format = :cmip,
    constraint_ts = (start_year = 2015, end_year = 2017) # for filename_format=:cmip, constraint_ts can be applied already for preview
)

# ---------------------------------- Load Model Data ------------------------------------- #
data_hist_proj_members = mwd.defineDataMap(
    paths_data, 
    data_ids; 
    level = :model,
    filename_format = :esmvaltool_cmip6,
    constraint = Dict(:models => ["AWI-CM-1-1-MR"]),
    constraint_ts = (start_year = 1950, end_year = 2014)
)



# ----------------------------- Load Observatioanl Data ---------------------------------- #
base_dir = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/obs/recipe_ERA5_20250718_180812/preproc/historical"
data_dirs = ["psl_CLIM-ann", "tas_CLIM-ann"]

obs_era5 = mwd.previewDataMap(
    joinpath.(base_dir, data_dirs),
    data_dirs;
    dtype = "observations",
    filename_format = :esmvaltool_obs,
    constraint = Dict(
        :variables => ["tas", "psl"]
    )
)

obs_era5 = mw.defineDataMap(
    joinpath.(base_dir, data_dirs),
    data_dirs;
    dtype = "observations",
    filename_format = :esmvaltool_obs,
    constraint = Dict(
        :variables => ["tas"]#, "psl"]
    ), 
    #constraint_ts = (start_year = 1950, end_year = 1991)
)

# ProfileSVG.@profview
# @profile
# @btime
#Profile.clear()
data_hist_proj_members = mwd.defineDataMap(
    paths_data, 
    data_ids; 
    level = :model,
    filename_format = :esmvaltool_cmip6,
    constraint = Dict(:model => ["AWI-CM-1-1-MR"]),
    #constraint_ts = (start_year = 1950, end_year = 2014)
)

df = mwd.loadPreprocData(meta_data; constraint_ts, dtype = "cmip")
ProfileSVG.@profview mwd.loadPreprocData(meta_data; constraint_ts, dtype = "cmip")


@btime data_hist_proj_members = mwd.defineDataMap(
    paths_data, 
    ["tas_annual_historical"], #, "tas_annual_ssp585"],# "psl_annual_historical"]; 
    level = :model,
    filename_format = :esmvaltool_cmip6,
    #constraint = Dict(:model => ["AWI-CM-1-1-MR"])
)



