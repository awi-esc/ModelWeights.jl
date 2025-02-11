import ModelWeights as mw
using NCDatasets
using DimensionalData
using Setfield


# LGM simulations
path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-anomalies-lgm-piControl.yml";
data =  mw.loadDataFromYAML(path_config, subset=Dict("level_shared_models" => mw.MEMBER));
reference = data["tas_CLIM_piControl"]
lgm_data = data["tas_CLIM_lgm"]
anomalies = mw.compute_anomalies(lgm_data, reference)

f1 = mw.makeSubplots(anomalies, (-20, 0), (nrows=3, ncols=5); figsize= (800,600).*(5,3))
save("plots/anomalies/lgm_anomalies.png", f1)

# Historical simulations
path_config = "/albedo/home/brgrus001/ModelWeights/configs/examples/example-anomalies-historical.yml";
data_all =  mw.loadDataFromYAML(path_config, subset=Dict("level_shared_models" => mw.MEMBER));
anomalies = mw.compute_anomalies(data_all["tas_CLIM_historical3"], data_all["tas_CLIM_historical0"])

# make plots for all members of a single model
members = Array(dims(anomalies, :member));
indices = findall(x -> startswith(x, "ACCESS-CM2#"), members)
data = anomalies[member = Where(x -> x in members[indices])]
joint_limits = (0, 2.5)
f2 = mw.makeSubplots(data, joint_limits, (nrows=2, ncols=5); figsize= (800,600).*(5,3))
save("plots/anomalies/historical_anomalies.png", f2)
