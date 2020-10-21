import JSON

function parse_truss_json(file_path::String)
    data = JSON.parsefile(file_path)
    ndim = data["dimension"]
    n = data["node_num"]
    m = data["element_num"]
    iT = Int
    T = Float64
    node_points = Dict{iT, SVector{ndim, T}}()
    elements = Dict{iT, Tuple{iT,iT}}()
    for i=1:n
        node_points[i] = data["node_list"][i]["point"]
        if "node_id" in keys(data["node_list"][i])
            @assert data["node_list"][i]["node_id"] == i
        end
    end
    for i=1:m
        elements[i] = (data["element_list"][i]...,)
    end
    E = data["E"]
    crosssecs = data["crosssecs"]
    if E isa Vector
        E = convert(Vector{T}, E)
        @assert length(E) == m
    else
        E = T(E)
    end
    if crosssecs isa Vector
        crosssecs = convert(Vector{T}, crosssecs)
        @assert length(crosssecs) == m
    else
        crosssecs = T(crosssecs)
    end
    return ndim, n, m, node_points, elements, E, crosssecs
end

function parse_support_load_json(file_path::String)
    data = JSON.parsefile(file_path)
    ndim = data["dimension"]
    iT = Int
    T = Float64

    n_load_nodes = length(data["point_load_list"])
    @assert(n_load_nodes > 0)
    loads = Dict{iT, SVector{ndim, T}}()
    for i=1:n_load_nodes
        load_v = data["point_load_list"][i]["node_id"]
        loads[load_v] = data["point_load_list"][i]["force"]
    end

    n_supp_nodes = length(data["support_node_list"])
    @assert(n_supp_nodes > 0)
    boundary = Dict{iT, SVector{ndim, Bool}}()
    for i=1:n_supp_nodes
        supp_v = data["support_node_list"][i]["node_id"]
        boundary[supp_v] = data["support_node_list"][i]["fixities"]
    end

    # TODO: self-weight
    return loads, boundary
end