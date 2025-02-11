import ModelWeights as mw
using DimensionalData
using Statistics
using Setfield
using NCDatasets
using CairoMakie


# Get data from piControl + historical + lgm experiments
path_config = "/albedo/home/brgrus001/ModelWeights/configs/lgm-historical.yml";
data_meta =  mw.loadDataFromYAML(path_config; preview = true, subset=Dict("level_shared_models" => mw.MODEL));
data = mw.loadDataFromYAML(path_config)

function makePlots(data, experiment, clim_var)
    members = dims(data.data, :member)
    for m in members
        f = Figure();
        mw.plotMeansOnMap!(f, data.data[member=At(m)], "$experiment mean $m")
        save("/albedo/home/brgrus001/ModelWeights/plots/$experiment-$clim_var/$m.png", f)
    end
end

makePlots(data["tos_CLIM_lgm"], "lgm", "tos")
makePlots(data["tos_CLIM_historical"], "historical", "tos")

