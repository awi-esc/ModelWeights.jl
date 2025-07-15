import ModelWeights as mw
import ModelWeights.Data as mwd
import ModelWeights.Timeseries as mwt
import ModelWeights.Weights as mww

using DimensionalData
using CairoMakie
using YAXArrays

dtype = "cmip"

path_config = "./configs/projection-plots.yml";
meta_data = mw.defineDataMap(path_config; preview=true)
data_all = mw.defineDataMap(path_config; dtype)
data_all =  mw.defineDataMap(path_config; dtype, constraint=Dict("level_shared" => "model"))

data = mw.summarizeMembers(data_all)
mw.summarizeMembers!(data)
mwd.apply!(data, mwt.filterTimeseries, 2015, 2100)


quantiles =  [0.167, 0.833]
function getDataProjectionPlot(dat::YAXArray, quantiles::Vector{<:Number})
    gms = mwd.globalMeans(dat)
    unweighted_avg = mww.weightedAvg(gms)
    uncertainties = mwd.uncertaintyRanges(gms; quantiles)
    return (avg=unweighted_avg, uncertainties=uncertainties)
end

historical = getDataProjectionPlot(data_ts["tas_CLIM-ann_historical"], quantiles);
ssp126 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp126"], quantiles);
ssp585 = getDataProjectionPlot(data_ts["tas_CLIM-ann_ssp585"], quantiles);

min_val = minimum(map(x -> minimum(x), [ssp126.avg, ssp585.avg]))
max_val = minimum(map(x -> maximum(x), [ssp126.avg, ssp585.avg]))


f = Figure(); 
ax = Axis(f[1,1], title = "Near-Surface Air Temperature", xlabel="Year");
ylims!(ax, min_val - 0.5, max_val + 5)

label_unc_unw(quantiles) = "Non-weighted quantiles: " * join(quantiles, "-")
label_unc_w(quantiles) = "Weighted quantiles: " * join(quantiles, "-")

mw.Plots.plotTimeseries!(ax, historical.avg; 
    uncertainties = coalesce.(historical.uncertainties, missing => NaN), 
    label="Non-weighted mean", 
    label_unc = label_unc_unw(historical.uncertainties), 
    color=:grey
)
mw.Plots.plotTimeseries!(
    ax, 
    ssp126.avg; 
    color_line=:blue, 
    uncertainties = coalesce.(ssp126.uncertainties, missing => NaN),
    label="Non-weighted mean", 
    label_unc = label_unc_unw(quantiles)
)
mw.Plots.plotTimeseries!(
    ax, 
    ssp585.avg; 
    uncertainties = coalesce.(ssp585.uncertainties, missing => NaN), 
    label = "Non-weighted mean", 
    label_unc = label_unc_unw(quantiles)
)

axislegend(position=:lt)
f