import JSON

function parse_truss_json(file_path::String)
    data = JSON.parsefile(file_path)
    ndim = data["dimension"]
    n = data["node_num"]
    m = data["element_num"]

    node_points = Dict{Int, SVector{ndim, Float64}}()
    elements = Dict{Int, Tuple{Int,Int}}()

    # get node coord
    for i=1:n
        node_points[i] = data["node_list"][i]["point"]
    end
    # get element node ids
    for i=1:m
        elements[i] = (data["element_list"][i]...,)
    end

    return ndim, n, m, node_points, elements
end

function parse_support_load_json(file_path::String)
    data = JSON.parsefile(file_path)
    ndim = data["dimension"]
    n = data["node_num"]

    n_load_nodes = length(data["point_load_list"])
    @assert(n_load_nodes > 0)
    loads = Dict{Int, SVector{ndim, Float64}}()
    for i=1:n_load_nodes
        load_v = data["point_load_list"][i]["node_id"]
        loads[load_v] = data["point_load_list"][i]["force"]
    end

    n_supp_nodes = length(data["support_node_list"])
    @assert(n_supp_nodes > 0)
    boundary = Dict{Int, SVector{ndim, Bool}}()
    for i=1:n_supp_nodes
        supp_v = data["support_node_list"][i]["node_id"]
        boundary[supp_v] = data["support_node_list"][i]["fixities"]
    end

    # TODO: self-weight
    return loads, boundary
end