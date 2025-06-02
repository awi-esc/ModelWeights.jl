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
    addDiagnostic!(datamap::DataMap, fn::Function, args...; kwargs..., id_suffix::String="")

Compute a diagnostic by running `fn` with positional arguments `args` and 
keyword arguments `kwargs` and add result to `datamap` at the id of the 
computed result concatenated with `#id_suffix` if `id_suffix` is not empty.
"""
function addDiagnostic!(
    datamap::DataMap, fn::Function, id::String, key_suffix::String, args...; kwargs...
)
    getAvailableIds(datamap, [id])
    data = fn(datamap[id], args...; kwargs...)
    new_key = isempty(key_suffix) ? data.properties["_id"] : join([data.properties["_id"], key_suffix], "#")
    if new_key in keys(datamap)
        @warn "$new_key overwritten!"
    end
    datamap[new_key] = data
    @info "run $(String(Symbol(fn))) for $id, added new entry at $new_key."
    return new_key
end


function addDiagnostic!(
    datamap::DataMap, fn::Function, ids::Vector{String}, key_suffix::String, args...; kwargs...
)
    ids = getAvailableIds(datamap, ids)
    new_keys = similar(ids)
    for (i, id) in enumerate(ids)
        new_keys[i] = addDiagnostic!(datamap, fn, id, key_suffix, args...; kwargs...)
    end
    return new_keys
end



"""
    addDiagnostic!(
        data::ClimateData, fn::Function, ids::Vector{String}, key_suffix::String, args...; kwargs...
    )
Run `addDiagnostic!` for model and observational data in `data`.
"""
function addDiagnostic!(
    data::ClimateData, fn::Function, ids::Vector{String}, key_suffix::String, args...; kwargs...
)
    @info "run addDiagnostic! for Model data..."
    addDiagnostic!(data.models, fn, ids, key_suffix, args...; kwargs...)
    @info "run addDiagnostic! for Observational data..."
    addDiagnostic!(data.obs, fn, ids, key_suffix, args...; kwargs...)
    return nothing
end