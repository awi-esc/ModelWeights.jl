import ModelWeights as mw
import ModelWeights.Data as mwd
import ModelWeights.Plots as mwp
using CairoMakie
using ColorSchemes
using DimensionalData
using Distributions
using LinearAlgebra
using Turing
using MCMCChains
using StatsPlots
using YAXArrays

plot_dir = "/albedo/home/brgrus001/ModelWeightsPaper/work/output/plots/likelihoods"
target_data_dir =  "/albedo/home/brgrus001/ModelWeightsPaper/work/output/data"
data = mwd.readDataFromDisk(joinpath(target_data_dir, "example-brunner_data.jld2"))

model_data = data["diagnostic_data"]["models"]
mwd.summarizeMembers!(model_data)
obs_data = data["diagnostic_data"]["observations"]

diagnostics = sort(["tas_ANOM", "tas_STD", "psl_ANOM", "psl_STD", "tas_TREND"])

# Plot Data
figs_rmse = []
figs_bias = []
for d in diagnostics
    @info "running plots for $d"
    distances = model_data[d] .- obs_data[d][model=1]
    s = size(distances)
    # colors = get(ColorSchemes.viridis, range(0,1; length=s[1]))
    # f = Figure()
    # ax = Axis(f[1,1], xlabel="Model bias $d")
    # for row in 1:s[1]
    #     for col in 1:s[2]
    #         Makie.density!(ax, distances[row,col,:], color=colors[row])
    #     end
    # end
    # push!(figs_bias, f)
    rmses = mwd.distancesData(model_data[d], obs_data[d]);
    mean_rmse = mean(rmses)
    xs = range(-2*mean_rmse, 2*mean_rmse; length=500)
    ys = Distributions.pdf.(Normal(0, mean_rmse), xs)
    f_rmse = Figure()
    ax = Axis(f_rmse[1,1], xlabel="RMSE $d")
    Makie.density!(ax, Array(rmses.data))
    Makie.scatter!(ax, xs, ys)
    push!(figs_rmse, f_rmse)
end

# use softmax weights based on RMSEs
weights = []
nb_models_above_avg = Dict()
for d in diagnostics
    rmses = mwd.distancesData(model_data[d], obs_data[d])
    models = collect(lookup(rmses, :model))
    w_unnormalized = exp.(-1 .* rmses)
    w = w_unnormalized ./ sum(w_unnormalized)
    w_yax = YAXArray((Dim{:model}(models), Dim{:weight}([d])), reshape(w, length(w), 1))
    push!(weights, w_yax)
    nb_models_above_avg[d] = sum(w .>= 1 / length(models))
end
weights = cat(weights...; dims=Dim{:weight}(diagnostics))
nb_models_above_avg

f_softmax = []
for d in diagnostics
    f = mwp.plotWeights(weights[weight = Where(x -> x == d)])
    push!(f_softmax, f)
end

climwip_weights = mwd.readDataFromDisk(joinpath(target_data_dir, "example-brunner_climwip-weights.jld2"))
f = Figure()
ax = Axis(f[1,1], xlabel="Genarlized distance Di")
Makie.density!(ax, vec(climwip_weights.Di.data))
f

f = Figure()
ax = Axis(f[1,1], xlabel="Normalized distances")
for d in diagnostics
    df = climwip_weights.performance_distances[d]
    df = mwd.summarizeMembers(df)
    df_normalized = df ./ median(df)
    Makie.density!(ax, collect(df_normalized), label="$d")
end
axislegend()
f



config = climwip_weights.config.performance
diagnostic_keys = collect(keys(config))
indices = sortperm(diagnostic_keys)
diagnostic_weights = YAXArray((Dim{:weight}(diagnostic_keys[indices]),), collect(values(config))[indices])

w_weighted_unnormalized = @d diagnostic_weights .* weights
w_weighted_unnormalized = sum(w_weighted_unnormalized; dims=:weight)[weight = 1]
w_weighted = w_weighted_unnormalized ./ sum(w_weighted)
w_weighted = YAXArray((Dim{:model}(collect(lookup(w_weighted, :model))), Dim{:weight}(["likelihood based weighted"])), reshape(vec(w_weighted), length(vec(w_weighted)), 1))

w_unnormalized = reduce(*, weights, dims=:weight)[weight=1]
w = w_unnormalized ./ sum(w_unnormalized)
w = YAXArray((Dim{:model}(collect(lookup(w, :model))), Dim{:weight}(["likelihood based"])), reshape(vec(w), length(vec(w)), 1))


all_weights = cat(climwip_weights.w[weight = Where(x -> x in [ "wP-historical"])], w_weighted, w; dims=:weight)
f = mwp.plotWeights(all_weights[weight = Where(x -> x in ["wP-historical", "likelihood based"])])
f = mwp.plotWeights(all_weights[weight = Where(x -> x in ["wP-historical", "likelihood based weighted"])])

# 1.use guassian error model for generalized distances, using different options for sigma
# the greater sigma, the more the weight is distributed equally across models
# the smaller sigma, the more the weight is centered on few models only
Di = climwip_weights.Di
indices = sortperm(Di)
Di = Di[indices]
models_sorted = collect(lookup(Di, :model))
indices_alphabet = sortperm(models_sorted)
generalized_distances = vec(Di.data)
f = Figure()
ax = Axis(f[1,1], xlabel="Genarlized distance Di")
for sigma_D in 0.02:0.1:1.92
    w_unnormalized = Distributions.pdf.(Normal(0, sigma_D/sqrt(2)), vec(Di.data))
    w = w_unnormalized ./ sum(w_unnormalized)
    Makie.scatterlines!(ax, generalized_distances, w, label="sigmaD=$sigma_D")
end
axislegend()
f

sigma_D = climwip_weights.config.sigma_performance
distr = Normal(0, sigma_D/sqrt(2))
w_unnormalized = map(x -> Distributions.pdf(distr, x), vec(Di.data))
w = w_unnormalized ./ sum(w_unnormalized)
w_normal = YAXArray((Dim{:model}(models_sorted[indices_alphabet]), Dim{:weight}(["Normal"])), reshape(w[indices_alphabet], length(w), 1))


# 2.use a multivariate normal distribution across diagnostics
data = Dict{String, YAXArray}(climwip_weights.performance_distances)
mwd.summarizeMembers!(data)
mwd.apply!(data, x -> x./median(x)) # normalize distances for each diagnostic

n_diagnostics = length(data)
n_models = length(data[diagnostics[1]])
data_mat = YAXArray((Dim{:model}(models_sorted), Dim{:diagnostic}(diagnostics)), rand(n_models, n_diagnostics))
for d in diagnostics
    data_mat[diagnostic = At(d)] .= vec(data[d][indices].data)
end

corr_mat = cor(data_mat; dims=:model)
corr_mat = Matrix(Symmetric(corr_mat))
epsilon = 1*10^-8;

Sigma = cov(data_mat; dims=:model)

f = Figure()
ax = Axis(f[1,1], xlabel="Genarlized distance Di")
for sigma_D in 0.02:0.1:1.92
    Sigma = (sigma_D/sqrt(2))^2 * (corr_mat + epsilon * I(n_diagnostics))
    distr = MvNormal(zeros(n_diagnostics), Sigma)
    w_unnormalized = map(x -> Distributions.pdf(distr, x), eachrow(Array(data_mat.data)))
    w = w_unnormalized ./ sum(w_unnormalized)
    Makie.scatterlines!(ax, generalized_distances, w, label="sigmaD=$sigma_D")
end
axislegend()
f

sigma_D = climwip_weights.config.sigma_performance
Sigma = (sigma_D)^2 * (corr_mat + epsilon * I(n_diagnostics))
distr = MvNormal(zeros(n_diagnostics), Sigma)
w_unnormalized = map(x -> Distributions.pdf(distr, x), eachrow(Array(data_mat.data)))
w = w_unnormalized ./ sum(w_unnormalized)
w_mvn = YAXArray((Dim{:model}(models_sorted[indices_alphabet]), Dim{:weight}(["MvNormal"])), reshape(w[indices_alphabet], length(w), 1))

weights_names = ["wP-historical", "Normal"]
weights_names = ["wIP-historical", "MvNormal"]
all_weights = cat(climwip_weights.w[weight = Where(x -> x in weights_names)], w_mvn, w_normal; dims=:weight)
f = mwp.plotWeights(all_weights)


# 






weights = []
nb_models_above_avg = Dict()
for d in diagnostics
    rmses = mwd.distancesData(model_data[d], obs_data[d])
    models = collect(lookup(rmses, :model))
    w_unnormalized = Distributions.pdf.(Normal(0, sigma_d/sqrt(2)), vec(rmses.data))
    w = w_unnormalized ./ sum(w_unnormalized)
    w_yax = YAXArray((Dim{:model}(models), Dim{:weight}([d])), reshape(w, length(w), 1))
    push!(weights, w_yax)
    nb_models_above_avg[d] = sum(w .>= 1 / length(models))
end
weights = cat(weights...; dims=Dim{:weight}(diagnostics))
nb_models_above_avg

all_weights = cat(climwip_weights.w[weight = Where(x -> x in [ "wP-historical"])], weights; dims=:weight)

f = mwp.plotWeights(all_weights[weight = Where(x -> x in ["wP-historical", diagnostics...])])




@model function likelihoodData()
    sigma ~ truncated(Normal(0, 1), 0, Inf)
    data ~ Normal(0, sigma)
end

MAPs_sigma = Dict()
for d in diagnostics
    @info "$d"
    distances = model_data[d] .- obs_data[d][model=1]
    data = vec(Array(distances.data))
    model = likelihoodData() | (data=data,);
    map_sigma = maximum_a_posteriori(model)
    sigma_hat = map_sigma.values[:sigma]
    MAPs_sigma[d] = sigma_hat
end





figures_data = []
fitted_distributions = Dict()
rmse_data = Dict()
distances = Dict()
for d in diagnostics
    distance = model_data[d] .- obs_data[d][model=1]      
    all_data = vec(Array(distance.data))
    # rmses = mwd.distancesData(model_data[d], obs_data[d]);
    # rmse_data[d] = rmses
    # all_data = vec(rmses.data)
    fitted = Distributions.fit_mle(Distributions.Normal, all_data)
    mu = fitted.μ
    sigma = fitted.σ
    
    
    f = Figure()
    ax = Axis(f[1,1], xlabel="$d", ylabel="Density")
    #Makie.density!(ax, all_data)
    #distance = reshape(Array(distance.data), :, size(distance)[3])
    Makie.density!(ax, all_data, color=:green, label="bias")
    
    
    xs_start, xs_end = extrema(all_data)
    xs = collect(range(xs_start, xs_end; length=1000))
    
    
    fitted_distr = Distributions.Normal(mu, sigma)
    Makie.scatter!(ax, xs, Distributions.pdf.(fitted_distr, xs), color=:steelblue)
    axislegend()
    f
    push!(figures_data, f)
    fitted_distributions[d] = fitted_distr
end


log_likelihoods = Dict()
weights = []
nb_models_above_avg = Dict()
for d in diagnostics
    data = rmse_data[d]
    models = collect(lookup(data, :model))
    likelihoods = Distributions.pdf.(fitted_distributions[d], vec(data.data))
    #log_likelihoods[d] = Vector(undef, length(models))
    # for (i, m) in enumerate(models)
    #     ll = sum(log.(Distributions.pdf.(fitted_distributions[d], data[model = At(m)].data[])))
    #     log_likelihoods[d][i] = ll
    # end
    #w = log_likelihoods[d]./sum(log_likelihoods[d])
    w = likelihoods ./ sum(likelihoods)
    w_yax = YAXArray((Dim{:model}(models), Dim{:weight}([d])), reshape(w, length(w), 1))
    push!(weights, w_yax)
    nb_models_above_avg[d] = sum(w .>= 1 / length(models))
end
weights = cat(weights...; dims=Dim{:weight}(diagnostics))
nb_models_above_avg

f = mwp.plotWeights(weights)



combined_w = YAXArray((Dim{:model}(collect(lookup(w, :model))), Dim{:weight}(["combined"])), reshape(vec(w), length(vec(w)), 1))



mwp.plotWeights(combined_w)
climwip_weights = mwd.readDataFromDisk(joinpath(target_data_dir, "example-brunner_climwip-weights.jld2"))

all_weights = cat(climwip_weights.w, combined_w; dims=:weight)
mwp.plotWeights(all_weights[weight = Where(x -> x in ["wP-historical", "combined"])])



MAPs_sigma = Dict()
# m = ones(2,3,4)
# m[:,2,:] .= 2
# m[:,3,:] .= 3
# reshape(m,:,4)
# mr = reshape(m, 6, 4)

@model function likelihoodData()
    sigma ~ truncated(Normal(0, 1), 0, Inf)
    data ~ Normal(0, sigma)
end

figures = []
for d in diagnostics    
    distance = model_data[d] .- obs_data[d][model=At("model1")]
    distance = mwd.summarizeMembers(distance)
    all_data = vec(Array((distance.data)))

    model = likelihoodData() | (data=all_data,);
    map_sigma = maximum_a_posteriori(model)
    
    models = lookup(distance, :model)
    palette = get(ColorSchemes.viridis, range(0,1; length=length(models)))
    f = Figure()
    ax = Axis(f[1,1], xlabel="$d: Model Bias", ylabel="Density")
    for m in models
        diffs = distance[model = At(m)]
        density!(ax, vec(diffs.data), alpha=0.8)
    end
    density!(ax, all_data)

    sigma_hat = map_sigma.values[:sigma]
    MAPs_sigma[d] = sigma_hat
    fitted_distr = Distributions.Normal(0.0, sigma_hat)
    xs_start, xs_end = extrema(all_data)
    xs = collect(range(xs_start, xs_end; length=1000))
    Makie.scatter!(ax, xs, Distributions.pdf.(fitted_distr, xs))
    push!(figures, f)
end


# Check RMSEs
MAPs_sigma_rmse = Dict()
figures = []
for d in diagnostics    
    rmses = mwd.distancesData(model_data[d], obs_data[d]);
    all_data = vec(rmses.data)

    # model = likelihoodData() | (data=all_data,);
    # map_sigma = maximum_a_posteriori(model)
    
    f = Figure()
    ax = Axis(f[1,1], xlabel="RMSE $d", ylabel="Density")
    Makie.density!(ax, all_data)

    # sigma_hat = map_sigma.values[:sigma]
    # MAPs_sigma_rmse[d] = sigma_hat
    # fitted_distr = Distributions.Normal(0.0, sigma_hat)
    # xs_start, xs_end = extrema(all_data)
    # xs = collect(range(xs_start, xs_end; length=1000))
    #Makie.scatter!(ax, xs, Distributions.pdf.(fitted_distr, xs))
    #Makie.density!(ax, rand(fitted_distr, 1000))
    push!(figures, f)
end

# compute model weights
log_likelihoods = Dict()
weights = []
nb_models_above_avg = Dict()
for d in diagnostics    
    likelihood_fn = Distributions.Normal(0, MAPs_sigma[d])
    distance = model_data[d] .- obs_data[d][model=At("model1")]
    distance = mwd.summarizeMembers(distance)
    models = collect(lookup(distance, :model))
    log_likelihoods[d] = Vector(undef, length(models))
    for (i,m) in enumerate(models)
        ll = sum(log.(Distributions.pdf.(likelihood_fn, vec(Array(distance[model = At(m)])))))
        log_likelihoods[d][i] = ll
    end
    w = log_likelihoods[d]./sum(log_likelihoods[d])
    w_yax = YAXArray((Dim{:model}(models), Dim{:weight}([d])), reshape(w, length(w), 1))
    push!(weights, w_yax)
    nb_models_above_avg[d] = sum(w .>= 1 / length(models))
end
weights = cat(weights...; dims=Dim{:weight}(diagnostics))

# plot weight vectors
f_weights = mwp.plotWeights(weights)




@model function weight(models)
    N_models = length(models)
    models ~ Categorical(fill(1/N_models, N_models))
    sigma ~ truncated(Normal(0, 1), 0, Inf)
    data ~ Normal(0, sigma)
end

for d in diagnostics    
    distance = model_data[d] .- obs_data[d][model=At("model1")]
    distance = mwd.summarizeMembers(distance)
    models = collect(lookup(distance, :model))

    all_data = vec(Array((distance.data)))
    model = weight(models) | (data=all_data,);
    chn = sample(model, NUTS(), 2_000, progress=false)

end



prior_samples = sample(distr, Prior(), 2_000, progress=false)
posterior_samples = sample(distr, NUTS(), 2_000, progress=false)

f_posterior = Figure()
ax = Axis(f_posterior[1,1], xlabel="sigma")
density!(ax, posterior_samples[:sigma].data[:,1])
f_posterior

distances = mwd.distancesData(model_data[diagnostic], obs_data[diagnostic]);
N_models = length(distances)

nb_models_above_avg = []
sigmas = collect(range(0.05, 2.05, 100))
for sigma in sigmas
    likelihoods = exp.(-(distances.data.^2 / (2 * sigma^2)))
    priors = 1/N_models
    weights_unnormalized = likelihoods * priors
    weights = weights_unnormalized/sum(weights_unnormalized)
    push!(nb_models_above_avg, sum(weights .>= 1/N_models))
end

distances = mwd.distancesData(model_data[diagnostic], obs_data[diagnostic]);
f = Figure()
ax = Axis(f[1,1], xlabel="RMSE $diagnostic")
Makie.density!(ax, distances.data)
f



diagnostics = map(x -> x * "-GM", ids)
mwd.apply!(model_data, mwd.globalMeans;
 ids = ids,
 ids_new = diagnostics
)
mwd.apply!(obs_data, mwd.globalMeans;
 ids = ids,
 ids_new = diagnostics
)


@model function weight(model_data, obs)
    N_models = length(model_data)
    sigma ~ truncated(Normal(0, 1), 0, Inf)
    #w ~ Dirichlet(fill(1.0, N_models))
    w = repeat([1 / N_models], N_models)
    #y_hat = sum(model_data .* w) # predicted weighted average
    obs ~ Normal(y_hat, sigma)

end

priors = Dict();
posteriors = Dict();
mean_estimates = Dict();
map_estimates = Dict();
for diagnostic in diagnostics
    @info "running $diagnostic ..."

    distances = mwd.distancesData(model_data[diagnostic], obs_data[diagnostic]);
    data = model_data[diagnostic]
    distr = weight(data.data,  obs_data[diagnostic].data[1]);
    chain_prior = sample(distr, Prior(), 2_000, progress=false)
    priors[diagnostic] = chain_prior
    chain_posterior = sample(distr,  NUTS(), 2_000, progress=false)
    posteriors[diagnostic] = chain_posterior

    f_posterior = plot(chain_posterior)
    mwp.savePlot(f_posterior, joinpath(plot_dir, diagnostic * "-posterior.png"))
    f_prior = plot(chain_prior)
    mwp.savePlot(f_prior, joinpath(plot_dir, diagnostic * "-prior.png"))

    sigma_samples = chain_posterior[:sigma].data  # extract sigma samples
    mean_estimates[diagnostic] = mean(sigma_samples)
    map_estimates[diagnostic] = maximum_a_posteriori(distr)

    N_models = length(data)
    w = repeat([1 / N_models], N_models)
    mu = sum(data.data .* w)
    likelihood = Distributions.Normal(mu, sigmas[diagnostic])
    f = Figure()
    ax = Axis(f[1,1], xlabel = "$diagnostic")
    lower_bound, upper_bound = mu - 5 * sigmas[diagnostic], mu + 5 * sigmas[diagnostic]
    xs = range(lower_bound, upper_bound; length = 1000)
    Makie.scatter!(ax, collect(xs), Distributions.pdf.(likelihood, xs))
    Makie.scatter!(ax, obs_data[diagnostic].data, [0]; color=:green, markersize=16, label="Observed")
    Makie.scatter!(ax, data.data, repeat([0], N_models); color=:orange, label="Models")
    Makie.scatter!(ax, [mean(data.data)], [0]; color=:indianred, markersize=16, label="Ensemble mean")
    axislegend()
    mwp.savePlot(f, joinpath(plot_dir, diagnostic * ".png"); overwrite = true)
end
#plot a normal distributionm

dists_perform_all = mwd.distancesData(model_data, obs_data, ["tas_TREND"]);

target_ecs_distr = Distributions.Gamma(67.696, 0.0476)




@model function ecsWeights(; Nmodels::Int)
    w ~ Dirichlet(fill(1.0, Nmodels))
    sqrt(sum(w .* models))
end

sampler = NUTS()
chain = sample(ecsWeights(;Nmodels=2), sampler, 2_000, progress=false)
histogram(chain)

ecsWeights(models) = ecsWeights(;Nmodels::Int = length(models)) | 




@model function gaussianWeight(model_data, obs)
    n_models = length(model_data)
    model ~ DiscreteUniform(1, n_models)
    sigma ~ truncated(Normal(0, 1), 0, Inf)
    obs ~ Normal(model_data[model], sigma^2)
end

sampler = NUTS();
chain = sample(gaussianWeight(
    model_data["tas_GM-ANOM"].data, 
    obs_data["tas_GM-ANOM"].data), sampler, 2_000, progress=false
)


@model function bmaWeights(models, obs)
    models = 
    Nmodels = length(models)
    w ~ Dirichlet(fill(1.0, Nmodels))
    sigma ~ truncated(Normal(0, 1), 0, Inf)

    yhat = sum(models .* w) # predicted weighted average
    obs ~ Normal(yhat, sigma^2)
    return w
end

bmaWeights(models, obs) = bmaWeights(models, obs) | (; obs)


sampler = NUTS()
prior = bmaWeights(model_data["tas_GM-ANOM"].data, obs_data["tas_GM-ANOM"].data)
chain = sample(prior, sampler, 2_000, progress=false)

histogram(chain)

# posterior 





# Unconditioned coinflip model with `N` observations.
@model function coinflip(; N::Int)
    # Our prior belief about the probability of heads in a coin toss.
    p ~ Beta(1, 1)

    # Heads or tails of a coin are drawn from `N` independent and identically
    # distributed Bernoulli distributions with success rate `p`.
    y ~ filldist(Bernoulli(p), N)

    return y
end;

data = coinflip(;N=10)
typeof(data)

rand(data)

coinflip(y::AbstractVector{<:Real}) = coinflip(; N=length(y)) | (; y)
obs = rand(data).y
model = coinflip(obs)

sampler = NUTS()

chain = sample(model, sampler, 2_000, progress=false);

histogram(chain)




@model function gaussian_mixture_model(x)
    # Draw the parameters for each of the K=2 clusters from a standard normal distribution.
    K = 2
    μ ~ MvNormal(Zeros(K), I)

    # Draw the weights for the K clusters from a Dirichlet distribution with parameters αₖ = 1.
    w ~ Dirichlet(K, 1.0)
    # Alternatively, one could use a fixed set of weights.
    # w = fill(1/K, K)

    # Construct categorical distribution of assignments.
    distribution_assignments = Categorical(w)

    # Construct multivariate normal distributions of each cluster.
    D, N = size(x)
    distribution_clusters = [MvNormal(Fill(μₖ, D), I) for μₖ in μ]

    # Draw assignments for each datum and generate it from the multivariate normal distribution.
    k = Vector{Int}(undef, N)
    for i in 1:N
        k[i] ~ distribution_assignments
        x[:, i] ~ distribution_clusters[k[i]]
    end

    return k
end
model = gaussian_mixture_model(x);
