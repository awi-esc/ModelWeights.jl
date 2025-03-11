import ModelWeights as mw
using DimensionalData
using CairoMakie
using YAXArrays

path_config = "./configs/projection-plots.yml";

meta_data = mw.loadDataFromYAML(path_config; preview=true)
data_all = mw.loadDataFromYAML(path_config)
# data_all =  mw.loadDataFromYAML(
#     path_config, subset=Dict("level_shared_models" => mw.MODEL)
# )
mw.averageEnsembleMembers!(data_all)
data_ts = mw.filterTimeseries(data_all, 2015, 2100)

# Brunner paper Fig. 2 with unweighted data
# required data: (Fig. 2a/2b for CMIP6/CMIP5)
# tas: ssp585 and ssp126 data until 2100
# tas: historical timeseries data

# compute ensemble means (unweighted and weighted)
# compute uncertainty ranges (weighted and unweighted)
# compute mean for historical (also weighted and unweighted)
function getDataProjectionPlot(dat::YAXArray)
    gms = mw.getGlobalMeansTS(dat)
    unweighted_avg = mw.computeWeightedAvg(gms)
    uncertainties = mw.getUncertaintyRanges(gms)
    return (avg=unweighted_avg, uncertainties=uncertainties)
end


historical = getDataProjectionPlot(data_ts["tas_CLIM-ann_historical"]);
ssp126 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp126"]);
ssp585 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp585"]);

f = Figure(); 
ax = Axis(f[1,1], title = "Near-Surface Air Temperature", xlabel="Year");

label_unc_unw(x) = "Non-weighted quantiles: " * join(x.properties["quantiles"], "-")
label_unc_w(x) = "Weighted quantiles: " * join(x.properties["quantiles"], "-")

mw.plotTimeseries!(f, ax, historical.avg; 
    uncertainties=coalesce.(historical.uncertainties, missing => NaN), 
    label="Non-weighted mean", label_unc=label_unc_unw(historical.uncertainties), 
    color=:grey)
mw.plotTimeseries!(ax, ssp126.avg; color=:blue, uncertainties=coalesce.(ssp126.uncertainties, missing => NaN),
    label="Non-weighted mean", label_unc=label_unc_unw(ssp126.uncertainties))
mw.plotTimeseries!(ax, ssp585.avg; uncertainties=coalesce.(ssp585.uncertainties, missing => NaN), 
    label="Non-weighted mean", label_unc=label_unc_unw(ssp585.uncertainties))
f

