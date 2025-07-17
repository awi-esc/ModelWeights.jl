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
data_all =  mw.defineDataMap(
    path_config; dtype, constraint=Dict("level_shared" => "model")
)
# or directly define start and end date:
data_all = mw.defineDataMap(
    path_config; 
    dtype, 
    constraint = Dict(
        "level_shared" => "model", 
        "timeseries" => Dict("start_y" => 2015, "end_y" => 2100)
    )
)

models = mwd.modelsFromMemberIDs(data_all["tas_CLIM-ann_ssp126"]; uniq = true)
#mwd.apply!(data_all, mwt.filterTimeseries, 2015, 2100; ids = ["tas_CLIM-ann_ssp126", "tas_CLIM-ann_ssp585"])
#mwd.apply!(data, mwt.filterTimeseries, -Inf, 2015; ids = ["tas_CLIM-ann_historical"])

data_historical = mw.defineDataMap(
    "/albedo/home/brgrus001/ModelWeights/configs/historical-annual.yml"; 
    dtype,
    constraint=Dict("level_shared" => "model", "models" => models, 
        "timeseries" => Dict("start_y" => 1850, "end_y" => 2014)
    )
)
models_historical =  mwd.modelsFromMemberIDs(data_historical["tas_CLIM-ann_historical"]; uniq = true)


quantiles =  [0.167, 0.833]
function getDataProjectionPlot(dat::YAXArray, quantiles::Vector{<:Number})
    gms = mwd.globalMeans(dat)
    unweighted_avg = mww.weightedAvg(gms)
    uncertainties = mwd.uncertaintyRanges(gms; quantiles)
    return (avg=unweighted_avg, uncertainties=uncertainties)
end

historical = getDataProjectionPlot(data_historical["tas_CLIM-ann_historical"], quantiles);
ssp126 = getDataProjectionPlot(data["tas_CLIM-ann_ssp126"], quantiles);
ssp585 = getDataProjectionPlot(data["tas_CLIM-ann_ssp585"], quantiles);

min_val = minimum(map(x -> minimum(x), [ssp126.avg, ssp585.avg, historical.avg]))
max_val = minimum(map(x -> maximum(x), [ssp126.avg, ssp585.avg, historical.avg]))


f = Figure(); 
ax = Axis(f[1,1], title = "Near-Surface Air Temperature", xlabel="Year");
ylims!(ax, min_val - 0.5, max_val + 5)

label_unc_unw(quantiles) = "Non-weighted quantiles: " * join(quantiles, "-")
label_unc_w(quantiles) = "Weighted quantiles: " * join(quantiles, "-")

mw.Plots.plotTimeseries!(
    ax, 
    historical.avg; 
    uncertainties = coalesce.(historical.uncertainties, missing => NaN), 
    label = "Non-weighted mean", 
    label_unc = label_unc_unw(quantiles), 
    color_line = :grey
)
mw.Plots.plotTimeseries!(
    ax, 
    ssp126.avg; 
    uncertainties = coalesce.(ssp126.uncertainties, missing => NaN),
    label = "Non-weighted mean", 
    label_unc = label_unc_unw(quantiles),
    color_line = :blue
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