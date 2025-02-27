import ModelWeights as mw
using DimensionalData
using Statistics
using Setfield
using CairoMakie
using ColorSchemes
using NCDatasets
using Distributions


# Get data from piControl + historical + lgm experiments
path_config = "/albedo/home/brgrus001/ModelWeights/configs/ecs-lgm-cooling.yml";
data_meta =  mw.loadDataFromYAML(path_config; preview = true);
data = mw.loadDataFromYAML(path_config; subset = Dict("level_shared_models" => mw.MODEL));

# for the shared models, make sure that physics of piControl models are the same
# as physics of lgm models
df = mw.alignPhysics(
    data, 
    data.map["tas_CLIM_lgm"].data.metadata["member_names"]; 
    level_shared_models = mw.MODEL
);
mw.averageEnsembleMembers!(df)
# lgm-cooling: here models need to be in same unit for both experiments
mw.kelvinToCelsius!(df)
lgm_cooling = df["tas_CLIM_lgm"].data .- df["tas_CLIM_piControl"].data;
# global lgm-cooling values for each model
# as we look at differences here, the unit (kelvin or celsius) doesnt matter
global_means = mw.getGlobalMeans(lgm_cooling)
predictions_gm_tas = Array(global_means);
ys = repeat([0], length(global_means));

# plot model data
begin
    f1 = Figure(size=(1000,400))
    t="GMST (LGM-piControl) with fitted normal distribution"
    ax = Axis(f1[1,1], title = L"$\Delta$ %$(t)")
    scatter!(ax, predictions_gm_tas, ys, color = :red)
    mu = Statistics.mean(global_means)
    sigma = Statistics.std(global_means; corrected=false, mean=mu)
    Distributions.fit_mle(Normal, global_means)

    dist = Distributions.Normal(mu, sigma)
    samples = rand(dist, 1000);
    hist!(samples)
    labels = Array(dims(global_means,:model))
    text!(ax, predictions_gm_tas, ys.+5, text=labels, rotation=pi/2)
    f1
end


# Assimilated data from Tierney et al. (2020)
tierney_data = NCDataset("/albedo/home/brgrus001/ModelWeightsPaper/work/data/Tierney-2020/Tierney2020_DA_atm.nc")
deltaSAT = DimArray(
    Array(tierney_data["deltaSAT"]), 
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"])))
)
errdeltaSAT = DimArray(
    Array(tierney_data["errdeltaSAT"]), 
    (Dim{:lon}(Array(tierney_data["lon"])), Dim{:lat}(Array(tierney_data["lat"])))
)
gm_deltaSAT = mw.getGlobalMeans(deltaSAT)
area_weights= mw.computeAreaWeights(Array(dims(deltaSAT, :lon)), Array(dims(deltaSAT, :lat)))
std_gm_deltaSAT = sum(area_weights .* errdeltaSAT)

CI_lower(n) = gm_deltaSAT - n * std_gm_deltaSAT;
CI_upper(n) = gm_deltaSAT + n * std_gm_deltaSAT; 

# 68% CI
l1 = CI_lower(1)
u1 = CI_upper(1)
# 95% CI
l2 = CI_lower(2)
u2 = CI_upper(2)

begin
    f3 = Figure(size = (1000, 400))
    t="Data assimilated reconstruction (Tierney et al., 2020)"
    ax = Axis(f3[1,1], title = L"%$(t) $\Delta$ GMST")
    xs = Array(global_means)
    ys = repeat([0],length(global_means))
    scatter!(ax, xs, ys)
    labels = Array(dims(global_means,:model))
    text!(ax, xs, ys.+5, text=labels, rotation=pi/2)
    dist = Distributions.Normal(gm_deltaSAT, std_gm_deltaSAT)
    samples = rand(dist, 1000);
    hist!(samples)
    vlines!(gm_deltaSAT; color=:grey, linewidth=5, label="Mean")
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
    f, ax = mw.makeScatterPlot(collect(1:length(global_means)), predictions_gm_tas;
        captions=(x="Models",y="Anomaly area-weighted global mean",title="LGM-cooling (lgm-piControl)"),
        xtick_labels = Array(dims(global_means, :model)),
        xticklabelrotation = pi/4,
        legend = (label="tas", position=:ct, color=:red),
        greyed_area = (
            y1=reconstruction.lower, y2=reconstruction.upper, 
            summary_val=reconstruction.summary_val, summary_stat=reconstruction.summary_stat, 
            label="$(reconstruction.ci) CI $(reconstruction.type) from Tierney et al. (2020)"
        )
    )
    save(reconstruction.target, f)
    f
end

###############################################################################
path_data = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/LGM";
path_recipes = "/albedo/home/brgrus001/ModelWeights/configs/lgm-cmip5-cmip6";
lgm_data = mw.loadDataFromESMValToolConfigs(
    path_data, path_recipes; 
    subset = Dict(
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
masksMissing = mapslices(x -> ismissing.(x), lgm_data.data, dims=:member);
longitudes = collect(dims(lgm_data.data, :lon));
latitudes = collect(dims(lgm_data.data, :lat));
members = collect(dims(lgm_data.data, :member));
area_weights = DimArray(zeros(length(longitudes), length(latitudes), length(members)), (Dim{:lon}(longitudes), Dim{:lat}(latitudes), Dim{:member}(members)));
units = lgm_data.data.metadata["units"];
for (idx, member) in enumerate(dims(lgm_data.data, :member))
    data = lgm_data.data[member = At(member)]
    mask = masksMissing[member = At(member)]
    area_weights[member = At(member)] = mw.computeAreaWeights(Array(dims(data, :lon)), Array(dims(data, :lat)); mask)
end
all(map(x -> area_weights[:,:, x] == area_weights[:,:,x+1], 1:length(members)-1))

# plot map of LGM tas data
fig1 = Figure()
mw.plotValsOnMap!(
    fig1,
    lgm_data.data[:,:,1], 
    "LGM tas for $(dims(lgm_data.data, :member)[1])";
    colors = ColorSchemes.twelvebitrainbow.colors   
)
fig2 = Figure();
mw.plotValsOnMap!(fig2, area_weights[:,:,1], "area weights for lgm models")


global_means_non_weighted = mapslices(x -> Statistics.mean(x), 
    lgm_data.data, dims=(:lon,:lat)
)[lon=1, lat=1]

global_means = mapslices(x -> Statistics.sum(x), 
    lgm_data.data .* area_weights, 
    dims=(:lon,:lat)
)[lon=1, lat=1]


obs_data = mw.loadDataFromESMValToolConfigs(
    "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_obs_historical_tas_tos", 
    "/albedo/home/brgrus001/ModelWeights/configs/historical_obs",
    dir_per_var = false,
    is_model_data = false,
    subset = Dict(
        "statistics" => ["CLIM"], 
        "variables" => ["tas"],
        "aliases" => ["historical"],
        "projects" => ["ERA5"] # for observational data default value is ["ERA5"]
    ),
    preview=false
)["tas_CLIM_historical"];


if obs_data.data.metadata["units"] == "K"
    obs_data = @set obs_data.data = mw.kelvinToCelsius(obs_data.data);
end
dat = obs_data.data[:,:,1]
fig3 = Figure();
mw.plotValsOnMap!(
    fig3, dat, 
    "Historical tas for $(dims(obs_data.data, :source)[1])",
    colors = ColorSchemes.twelvebitrainbow.colors
)
mask_obs = ismissing.(dat)
area_weights_obs = mw.computeAreaWeights(Array(dims(dat, :lon)), Array(dims(dat, :lat)); mask=mask_obs)
fig4 = Figure();
mw.plotValsOnMap!(fig4, area_weights_obs, "area weights for observations")
global_means_obs = mapslices(
    x -> Statistics.sum(skipmissing(x)), 
    dat .* area_weights_obs, 
    dims=(:lon,:lat)
)[lon=1, lat=1]


begin 
    f = Figure(size=(800,600)); 
    models = Array(dims(lgm_data.data, :member))
    ax = Axis(f[1,1], 
        xticks = (1:length(models), models), 
        #limits = ((years[1]-10, years[end]+10), (-1, 5)),
        title = "LGM global means",
        xlabel = "Model", 
        ylabel = "Gloabal Mean Temperature in °C",
        xticklabelrotation = pi/2
    );
    lines!(ax, 1:length(models), global_means.data, color = :red, label = "Area-weighted")
    scatter!(ax, 1:length(models), global_means.data, color = :red, label = "Area-weighted")

    lines!(ax, 1:length(models), global_means_non_weighted.data, color = :blue, label = "Non-weighted")
    scatter!(ax, 1:length(models), global_means_non_weighted.data, color = :blue, label = "Non-weighted")

    lines!(ax, 1:length(models), repeat([global_means_obs], length(models)), color = :green, label = "Area-weighted Global Mean Historical observations")
    axislegend(ax, merge = true, position = :cc)
    f
end


