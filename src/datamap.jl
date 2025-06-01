# Functions for adding data to DataMap objects.

function getAvailableIds(data::DataMap, ids::Vector{String})
    indices_keys = map(id -> haskey(data, id), ids)
    ids_absent = ids[.!(indices_keys)]
    if !isempty(ids_absent)
        if length(ids_absent) == length(ids)
            throw(ArgumentError("No data found for ids: $(ids)"))
        else
            @warn "No data found for ids" ids_absent
        end
    end
    return ids[indices_keys]
end



"""
    addDiagnostic!(datamap::DataMap, fn::Function, args...; kwargs...)

Compute a diagnostic by running `fn` with positional arguments `args` and 
keyword arguments `kwargs` and add result to `datamap` at the id of the 
computed result.
"""
function addDiagnostic!(datamap::DataMap, fn::Function, id::String, args...; kwargs...)   
    getAvailableIds(datamap, [id])
    data = fn(datamap[id], args...; kwargs...)
    @info "run $(String(Symbol(fn))) for $id, add new id $(data.properties["_id"])."
    datamap[data.properties["_id"]] = data
    return nothing
end


function addDiagnostic!(datamap::DataMap, fn::Function, ids::Vector{String}, args...; kwargs...)
    ids = getAvailableIds(datamap, ids)
    data_ids = collect(keys(datamap))
    new_ids = similar(ids)
    for (i, id) in enumerate(ids)
        new_id = fn(datamap[id], args...; kwargs..., only_newid=true)        
        if new_id in data_ids
            @warn "$new_id skipped computation since already present in datamap."
        else
            addDiagnostic!(datamap, fn, id, args...; kwargs...)
        end
        new_ids[i] = new_id
    end
    return new_ids
end


function addDiagnostic!(data::ClimateData, fn::Function, ids::Vector{String}, args...; kwargs...)
    @info "run addDiagnostic! for Model data..."
    addDiagnostic!(data.models, fn, ids, args...; kwargs...)
    @info "run addDiagnostic! for Observational data..."
    addDiagnostic!(data.obs, fn, ids, args...; kwargs...)
    return nothing
end