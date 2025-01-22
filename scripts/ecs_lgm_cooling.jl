import ModelWeights as mw
using DimensionalData
using Statistics
using Setfield
using CairoMakie
using ColorSchemes
using NCDatasets


# Get data from piControl + historical + lgm experiments
path_config = "/albedo/home/brgrus001/ModelWeights/configs/ecs-lgm-cooling.yml";
data_meta =  mw.loadDataFromYAML(path_config; preview = true);
data = mw.loadDataFromYAML(path_config; level_shared_models = mw.MODEL);

# for the shared models, make sure that physics of piControl models are the same
# as physics of lgm models
df = mw.alignPhysics(
    data, 
    data["tas_CLIM_lgm"].data.metadata["member_names"]; 
    level_shared_models = mw.MODEL
);
mw.averageEnsembleMembers!(df)

# lgm-cooling: here models need to be in same unit for both experiments
df_celsius = mw.kelvinToCelsius(df)

lgm_cooling = df_celsius["tas_CLIM_lgm"].data .- df_celsius["tas_CLIM_piControl"].data;

# global lgm-cooling values for each model
# as we look at differences here, the unit (kelvin or celsius) doesnt matter
global_means = mw.getGlobalMeans(lgm_cooling)

f1 = Figure()
ax = Axis(f1[1,1])
density!(global_means, direction=:x)
scatter!(ax, Array(global_means), repeat([0],length(global_means)), color = :red)
f1

f, ax = mw.makeScatterPlot(collect(1:length(global_means)), Array(global_means);
    captions=(x="Models",y="Global Mean Difference",title="LGM-cooling (lgm-piControl)"),
    xtick_labels = Array(dims(global_means, :model)),
    xticklabelrotation = pi/4,
    legend = (label="tas", position=:ct, color=:red),
    greyed_area = (
        y1=-6.8, y2=-4.4, 
        summary_val=-5.6, summary_stat="median", 
        label="95% confidence Proxy-only from Tierney et al. (2020)"
    )
)
f

# assimilated data from Tierney et al. (2020)
tierney_data = NCDataset("/albedo/home/brgrus001/ModelWeightsPaper/work/data/Tierney2020_DA_atm.nc")
deltaSAT = tierney_data["deltaSAT"]








path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; 
    subset = mw.Constraint(
        statistics = ["CLIM"], 
        variables = ["tos", "tas"],
        projects = ["CMIP5", "CMIP6"],
        aliases = ["lgm"], 
        subdirs = ["20241114"]
    )
)
lgm_tas = lgm_data["tas_CLIM_lgm"];
lgm_tos = lgm_data["tos_CLIM_lgm"];

lgm_data = lgm_tas;
# compute area weights for each model (depends on missing values)
masksMissing = mapslices(x -> ismissing.(x), lgm_data.data, dims=:member);
longitudes = collect(dims(lgm_data.data, :lon));
latitudes = collect(dims(lgm_data.data, :lat));
members = collect(dims(lgm_data.data, :member));
area_weights = DimArray(zeros(length(longitudes), length(latitudes), length(members)), (Dim{:lon}(longitudes), Dim{:lat}(latitudes), Dim{:member}(members)));
units = lgm_data.data.metadata["units"];
for (idx, member) in enumerate(dims(lgm_data.data, :member))
    data = lgm_data.data[member = At(member)]
    if (isa(units, Vector) && units[idx] == "K") || units == "K"
        df = lgm_data.data
        df[member = At(member)] = mw.kelvinToCelsius(df[member = At(member)])
        lgm_data = @set lgm_data.data = df
    end
    mask = masksMissing[member = At(member)]
    area_weights[member = At(member)] = mw.computeAreaWeights(data; mask)
end
all(map(x -> area_weights[:,:, x] == area_weights[:,:,x+1], 1:length(members)-1))

# plot map of LGM tas data
mw.plotMeansOnMap(
    lgm_data.data[:,:,1], 
    "LGM tas for $(dims(lgm_data.data, :member)[1])";
    colors = ColorSchemes.twelvebitrainbow.colors   
)
mw.plotMeansOnMap(area_weights[:,:,1], "area weights for lgm models")


global_means_non_weighted = mapslices(x -> Statistics.mean(x), 
    lgm_data.data, dims=(:lon,:lat)
)[lon=1, lat=1]


global_means = mapslices(x -> Statistics.sum(x), 
    lgm_data.data .* area_weights, 
    dims=(:lon,:lat)
)[lon=1, lat=1]


obs_data = mw.loadDataFromESMValToolConfigs(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_20241119_124434", 
    "/albedo/home/brgrus001/ModelWeights/configs/historical_obs",
    dir_per_var = false,
    is_model_data = false,
    subset = mw.Constraint(
        statistics = ["CLIM"], 
        variables = ["tas"],
        aliases = ["historical"],
        projects = ["ERA5"] # for observational data default value is ["ERA5"]
    )
)["tas_CLIM_historical"];


if obs_data.data.metadata["units"] == "K"
    obs_data = @set obs_data.data = mw.kelvinToCelsius(obs_data.data);
end
mw.plotMeansOnMap(
    obs_data.data[:,:,1], 
    "Historical tas for $(dims(obs_data.data, :source)[1])",
    colors = ColorSchemes.twelvebitrainbow.colors
)


mask_obs = ismissing.(obs_data.data[:,:,1])
area_weights_obs = mw.computeAreaWeights(obs_data.data[:,:,1]; mask)
mw.plotMeansOnMap(area_weights_obs, "area weights for observations")
global_means_obs = mapslices(
    x -> Statistics.sum(skipmissing(x)), 
    obs_data.data[:,:,1] .* area_weights_obs, 
    dims=(:lon,:lat)
)[lon=1, lat=1]




begin 
    f = Figure(); 
    models = Array(dims(lgm_data.data, :member))
    ax = Axis(f[1,1], 
        xticks = (1:length(models), models), 
        #limits = ((years[1]-10, years[end]+10), (-1, 5)),
        title = "LGM global means",
        xlabel = "Model", 
        ylabel = "Gloabal Mean Temperature in °C",
        xticklabelrotation = pi/2
    );
    lines!(ax, 1:length(models), vec(global_means), color = :red, label = "Area-weighted")
    scatter!(ax, 1:length(models), vec(global_means), color = :red, label = "Area-weighted")

    lines!(ax, 1:length(models), vec(global_means_non_weighted), color = :blue, label = "Non-weighted")
    scatter!(ax, 1:length(models), vec(global_means_non_weighted), color = :blue, label = "Non-weighted")

    scatter!(ax, 1, global_means_obs, color = :green, label = "Weighted Global Mean Observations")
    axislegend(ax, merge = true, position = :rc)
    f
end


