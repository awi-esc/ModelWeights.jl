import ModelWeights as mw
import ModelWeights.Data as mwd
import ModelWeights.Plots as mwp

using DimensionalData
using Statistics
using Setfield
using CairoMakie
using ColorSchemes
using NCDatasets
using Distributions
using YAXArrays
using Missings

# Get data from piControl + historical + lgm experiments
path_config = "./configs/ecs-lgm-cooling.yml";
data = mw.defineDataMap(
    path_config; 
    constraint = Dict("level_shared" => "model"),
    dtype = "cmip",
    filename_format = :esmvaltool
)


# for the shared models, make sure that physics of piControl models are the same
# as physics of lgm models
mwd.apply!(
    data, mwd.alignPhysics, Array(data["tas_CLIM_lgm"].member); 
    ids = ["tas_CLIM_piControl"],
    ids_new = ["tas_CLIM_piControl-lgm-members"]
)
mwd.summarizeMembers!(data)
mwd.apply!(
    data, mwd.anomalies, data["tas_CLIM_piControl"]; 
    ids = ["tas_CLIM_lgm"],
    ids_new = ["tas_ANOM_lgm-piControl"]
)





# NaN important for plotting! Type mustn't be Missing
global_means_tas = coalesce.(mwd.globalMeans(data["tas_ANOM_lgm-piControl"]), NaN);


# Assimilated data from Tierney et al. (2020)
global_means = global_means_tas;
tierney_data = NCDataset("/albedo/home/brgrus001/ModelWeightsPaper/work/data/Tierney-2020/Tierney2020_DA_atm.nc")
deltaSAT = YAXArray(
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"]))),
    Array(tierney_data["deltaSAT"])
)
errdeltaSAT = YAXArray(
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"]))),
    Array(tierney_data["errdeltaSAT"])
)
gm_delta = mwd.globalMeans(deltaSAT)[1]
area_weights_mat = mwd.makeAreaWeightMatrix(Array(dims(deltaSAT, :lon)), Array(dims(deltaSAT, :lat)))
std_gm_delta = sum(area_weights_mat .* errdeltaSAT)


# Only proxy data from Tierney et al. (2020)
global_means = global_means_tos;
tierney_data = NCDataset("/albedo/home/brgrus001/ModelWeightsPaper/work/data/Tierney-2020/Tierney2020_ProxyData_5x5_deltaSST.nc")
mat = allowmissing(Array(tierney_data["deltaSST"]))
mat[isnan.(mat)] .= missing
deltaSST = YAXArray(
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"]))),
    mat
)
mat = allowmissing(Array(tierney_data["std"]))
mat[isnan.(mat)] .= missing
errdeltaSST = YAXArray(
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"]))),
    mat
)
gm_delta = mw.globalMeans(deltaSST)[1]
# also adapt area weight matrix for missing values..!!?
area_weights_mat = mw.makeAreaWeightMatrix(
    Array(dims(deltaSST, :lon)), Array(dims(deltaSST, :lat));
    mask=ismissing.(mat)
)
std_gm_delta = sum(skipmissing(area_weights_mat .* errdeltaSST))




begin
    CI_lower(n) = gm_delta - n * std_gm_delta;
    CI_upper(n) = gm_delta + n * std_gm_delta; 

    # 68% CI
    l1 = CI_lower(1)
    u1 = CI_upper(1)
    # 95% CI
    l2 = CI_lower(2)
    u2 = CI_upper(2)
end




begin
    f3 = Figure(size = (1000, 400))
    t="Data assimilated reconstruction (Tierney et al., 2020)"
    ax = Axis(f3[1,1], title = L"%$(t) $\Delta$ GMST")
    xs = Array(global_means)
    ys = repeat([0],length(global_means))
    scatter!(ax, xs, ys)
    labels = Array(dims(global_means,:model))
    text!(ax, xs, ys.+5, text=labels, rotation=pi/2)
    dist = Distributions.Normal(gm_delta, std_gm_delta)
    samples = rand(dist, 1000);
    hist!(samples)
    vlines!(gm_delta; color=:grey, linewidth=5, label="Mean")
    lines!([l1, u1], [0, 0]; linewidth=5, color=:red, label="68% CI")
    lines!([l2, u2], [0, 0]; linewidth=5, color=:green, alpha=0.5, label="95% CI")
    axislegend(ax, merge = true, position = :rt)
    f3
end
save("plots/ecs-lgm-cooling/DA_GMST-lgm-cooling.png", f3)


# values from Tierney paper, DA-reproduced (CI-interval mismatch..!)
reconstruction = (type="proxy-only", lower=-6.8, upper=-4.4, summary_stat="median(?)", 
    summary_val=-5.6, ci="95%", target="plots/ecs-lgm-cooling/lgm_cooling_global_mean_DA.png");
reconstruction = (type="Data assimilated", lower=-6.5, upper=-5.7, summary_stat="mean", 
    summary_val=-6.1, ci="95%(68%??)", target="plots/ecs-lgm-cooling/lgm_cooling_global_mean_proxy_only.png");

begin
    f, ax = mwp.makeScatterPlot(collect(1:length(global_means_tas)), global_means_tas;
        captions=(x="Models",y="Anomaly area-weighted global mean",title="LGM-cooling (lgm-piControl)"),
        xtick_labels = Array(dims(global_means_tas, :model)),
        xticklabelrotation = pi/4,
        legend = (label="tas", position=:rc, color=:red),
        greyed_area = (
            y1=reconstruction.lower, y2=reconstruction.upper, 
            summary_val=reconstruction.summary_val, summary_stat=reconstruction.summary_stat, 
            label="$(reconstruction.ci) CI $(reconstruction.type) from Tierney et al. (2020)"
        )
    )
    #save(reconstruction.target, f)
    f
end

###############################################################################
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "./configs/lgm-cmip5-cmip6";
lgm_data = mw.loadDataFromESMValToolRecipes(
    path_data, path_recipes; 
    constraint = Dict(
        "statistics" => ["CLIM"], 
        "variables" => ["tos", "tas"],
        "projects" => ["CMIP5", "CMIP6"],
        "aliases" => ["lgm"], 
        "subdirs" => ["20241114"]
    )
)
mw.kelvinToCelsius!(lgm_data);
lgm_tas = lgm_data["tas_CLIM_lgm"];
lgm_tos = lgm_data["tos_CLIM_lgm"];

# just use tas data:
lgm_data = lgm_tas;
# compute area weights for each model (depends on missing values)
masksMissing = mapslices(x -> ismissing.(x), DimArray(lgm_data), dims=:member);
longitudes = collect(dims(lgm_data, :lon));
latitudes = collect(dims(lgm_data, :lat));
members = collect(dims(lgm_data, :member));
area_weights = DimArray(
    zeros(length(longitudes), length(latitudes), length(members)), 
    (Dim{:lon}(longitudes), Dim{:lat}(latitudes), Dim{:member}(members))
);
units = lgm_data.properties["units"];
for (idx, member) in enumerate(dims(lgm_data, :member))
    data = lgm_data[member = At(member)]
    mask = masksMissing[member = At(member)]
    area_weights[member = At(member)] = mw.makeAreaWeightMatrix(Array(dims(data, :lon)), Array(dims(data, :lat)); mask)
end
all(map(x -> area_weights[:,:, x] == area_weights[:,:,x+1], 1:length(members)-1))

# plot map of LGM tas data
fig1 = Figure();
mw.plotValsOnMap!(
    fig1,
    lgm_data[:,:,1], 
    "LGM tas for $(dims(lgm_data, :member)[1])";
    colors = ColorSchemes.twelvebitrainbow.colors   
)
fig1
fig2 = Figure();
mw.plotValsOnMap!(fig2, area_weights[:,:,1], "area weights for lgm models")
fig2

global_means_non_weighted = mapslices(x -> Statistics.mean(x), 
    DimArray(lgm_data), dims=(:lon,:lat)
)[lon=1, lat=1]

global_means = mapslices(x -> Statistics.sum(x), 
    DimArray(lgm_data) .* area_weights, 
    dims=(:lon,:lat)
)[lon=1, lat=1]


# TODO: config file non existent!
obs_data = mw.loadDataFromESMValToolRecipes(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_tas_tos", 
    "./configs/historical_obs",
    dir_per_var = false,
    is_model_data = false,
    constraint = Dict(
        "statistics" => ["CLIM"], 
        "variables" => ["tas"],
        "aliases" => ["historical"],
        "projects" => ["ERA5"] # for observational data default value is ["ERA5"]
    ),
    preview=false
)["tas_CLIM_historical"];


if obs_data.properties["units"] == "K"
    obs_data = mw.kelvinToCelsius(obs_data)
end
dat = obs_data[:,:,1]
fig3 = Figure();
mw.plotValsOnMap!(
    fig3, dat, 
    "Historical tas for $(dims(obs_data, :model)[1])",
    colors = ColorSchemes.twelvebitrainbow.colors
)
mask_obs = ismissing.(dat)
area_weights_obs = mw.makeAreaWeightMatrix(Array(dims(dat, :lon)), Array(dims(dat, :lat)); mask=mask_obs)
fig4 = Figure();
mw.plotValsOnMap!(fig4, area_weights_obs, "area weights for observations")
global_means_obs = mapslices(
    x -> Statistics.sum(skipmissing(x)), 
    dat .* area_weights_obs, 
    dims=(:lon,:lat)
)[lon=1, lat=1]


begin 
    f = Figure(size=(800,600)); 
    models = Array(dims(lgm_data, :member))
    ax = Axis(f[1,1], 
        xticks = (1:length(models), models), 
        #limits = ((years[1]-10, years[end]+10), (-1, 5)),
        title = "LGM global means",
        xlabel = "Model", 
        ylabel = "Gloabal Mean Temperature in °C",
        xticklabelrotation = pi/2
    );
    lines!(ax, 1:length(models), global_means, color = :red, label = "Area-weighted")
    scatter!(ax, 1:length(models), global_means, color = :red, label = "Area-weighted")

    lines!(ax, 1:length(models), global_means_non_weighted, color = :blue, label = "Non-weighted")
    scatter!(ax, 1:length(models), global_means_non_weighted, color = :blue, label = "Non-weighted")

    lines!(ax, 1:length(models), repeat([global_means_obs], length(models)), color = :green, label = "Area-weighted Global Mean Historical observations")
    axislegend(ax, merge = true, position = :cc)
    f
end


