using NCDatasets

climVar = "msftmz"
dim_key = "lev"

p = "/albedo/work/projects/p_forclima/preproc_data_esmvaltool/historical/recipe_historical_" * climVar * "/preproc/climatology_historical1/" * climVar;
files = readdir(p; join=true);
ncFiles = filter(x->endswith(x, ".nc"), files);
values = []
for file in ncFiles
    ds = NCDataset(file)
    dsVar = ds[climVar];
    push!(values, collect(dsVar[dim_key][:]));
end

identicals = unique(values)

indices1 = findall(x-> x==identicals[1], values)
indices2 = findall(x-> x==identicals[2], values)
indices3 = findall(x-> x==identicals[3], values)
indices4 = findall(x-> x==identicals[4], values)
indices13 = findall(x-> x==identicals[13], values)

filenames = readdir(p);
filenames= filter(x->endswith(x, ".nc"), filenames);
map(x -> split(x, ".")[1], filenames[indices2])