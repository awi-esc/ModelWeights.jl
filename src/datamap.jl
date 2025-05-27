# Functions for adding data to DataMap objects.

function errorIfIdsNotPresent(data::DataMap, ids::Vector{String})
    ismissing = map(x -> !haskey(data, x), ids)
    if any(ismissing)
        throw(ArgumentError("There is data missing for: $(ids[ismissing])!"))
    end
    return nothing
end


"""
    addDiagnostic!(datamap::DataMap, fn::Function, args...; kwargs...)

Compute a diagnostic by running `fn` with positional arguments `args` and 
keyword arguments `kwargs` and add result to `datamap` at the id of the 
computed result.
"""
function addDiagnostic!(datamap::DataMap, fn::Function, id::String, args...; kwargs...)   
    errorIfIdsNotPresent(datamap, [id])
    data = fn(datamap[id], args...; kwargs...)
    @info "run $(String(Symbol(fn))) for $id, add new id $(data.properties["_id"])."
    datamap[data.properties["_id"]] = data
    return nothing
end


function addDiagnostic!(datamap::DataMap, fn::Function, ids::Vector{String}, args...; kwargs...)
    errorIfIdsNotPresent(datamap, ids)
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
