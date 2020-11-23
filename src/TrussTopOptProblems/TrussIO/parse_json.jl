import JSON

function parse_truss_json(file_path::String)
    data = JSON.parsefile(file_path)
    ndim = data["dimension"]
    n = data["node_num"]
    m = data["element_num"]
    iT = Int
    T = Float64

    node_points = Dict{iT, SVector{ndim, T}}()
    for i=1:n
        node_points[i] = convert(SVector{ndim,T}, data["nodes"][i]["point"])
        if "node_ind" in keys(data["nodes"][i])
            @assert data["nodes"][i]["node_ind"] == i
        end
    end
    @assert length(node_points) == n

    elements = Dict{iT, Tuple{iT,iT}}()
    element_inds_from_tag = Dict()
    for i=1:m
        elements[i] = (data["elements"][i]["end_node_inds"]...,)
        if "elem_ind" in keys(data["elements"][i])
            @assert data["elements"][i]["elem_ind"] == i
        end
        elem_tag = data["elements"][i]["elem_tag"]
        if elem_tag ∉ keys(element_inds_from_tag)
            element_inds_from_tag[elem_tag] = []
        end
        push!(element_inds_from_tag[elem_tag], i)
    end
    @assert length(elements) == m

    E_from_tag = Dict()
    for mat in data["materials"]
        mat["elem_tags"] = length(mat["elem_tags"]) == 0 ? [nothing] : mat["elem_tags"]
        for e_tag in mat["elem_tags"]
            if e_tag in keys(E_from_tag)
                @warn "Multiple materials assigned to the same element tag |$(e_tag)|!"
            end
            E_from_tag[e_tag] = mat["E"]
        end
    end
    A_from_tag = Dict()
    for cs in data["cross_secs"]
        cs["elem_tags"] = length(cs["elem_tags"]) == 0 ? [nothing] : cs["elem_tags"]
        for e_tag in cs["elem_tags"]
            if e_tag in keys(A_from_tag)
                @warn "Multiple cross secs assigned to the same element tag |$(e_tag)|!"
            end
            A_from_tag[e_tag] = cs["A"]
        end
    end
    Es = zeros(T, m)
    As = zeros(T, m)
    for (tag, e_ids) in element_inds_from_tag
        for ei in e_ids
            if tag ∉ keys(E_from_tag)
                # use default material (key `nothing`)
                Es[ei] = E_from_tag[nothing]
                As[ei] = A_from_tag[nothing]
            end
        end
    end

    n_supp_nodes = length(data["supports"])
    # TODO only translation dof for now
    @assert(n_supp_nodes > 0)
    fixities = Dict{iT, SVector{ndim, Bool}}()
    for i=1:n_supp_nodes
        supp_v = iT(data["supports"][i]["node_ind"])
        fixities[supp_v] = data["supports"][i]["condition"]
    end

    load_cases = Dict()
    for (lc_ind, lc_data) in data["loadcases"]
        nploads = length(lc_data["ploads"])
        @assert nploads > 0
        ploads = Dict{iT, SVector{ndim, T}}()
        for pl in lc_data["ploads"]
            load_v = pl["node_ind"]
            ploads[load_v] = convert(SVector{ndim,T}, pl["force"])
        end
        load_cases[lc_ind] = ploads
    end

    return node_points, elements, Es, As, fixities, load_cases
end