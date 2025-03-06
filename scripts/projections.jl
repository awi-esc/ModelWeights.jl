import ModelWeights as mw
using DimensionalData
using CairoMakie

path_config = "configs/projection-data.yml";

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
function getDataProjectionPlot(dat::mw.Data)
    gms = mw.getGlobalMeansTS(dat.data)
    unweighted_avg = mw.computeWeightedAvg(gms)
    uncertainties = mw.getUncertaintyRanges(gms)

    unweighted_avg.properties["label"] = dat.meta.attrib.exp * ": Non-weighted mean"
    uncertainties.properties["label"] = dat.meta.attrib.exp * ": Non-weighted quantiles: " * join(uncertainties.properties["quantiles"], "-")    
    return (avg=unweighted_avg, uncertainties=uncertainties)
end


historical = getDataProjectionPlot(data_ts["tas_CLIM-ann_historical"]);
ssp126 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp126"]);
ssp585 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp585"]);

f = Figure(); 
ax = Axis(f[1,1], title = "Near-Surface Air Temperature", xlabel="Year");
mw.plotTimeseries!(f, ax, historical.avg; uncertainties = historical.uncertainties, color=:grey)
mw.plotTimeseries!(f, ax, ssp126.avg; color=:blue, uncertainties = ssp126.uncertainties)
mw.plotTimeseries!(f, ax, ssp585.avg; uncertainties = ssp585.uncertainties)
f



# Check pseudo observations in their Fig.2!

