using TopOpt.TopOptProblems: StiffnessTopOptProblem, Metadata

@params struct TrussProblem{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
    truss_grid::TrussGrid{dim, T, N, M}
    E::T
    ν::T
    ch::ConstraintHandler{<:DofHandler{dim, N, T, M}, T}
    black::AbstractVector
    white::AbstractVector
    varind::AbstractVector{Int}
    metadata::Metadata
end

function TrussProblem(::Type{Val{CellType}}, node_points::Dict{iT, SVector{dim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
    loads::Dict{iT, SVector{dim, T}}, supports::Dict{iT, SVector{dim, fT}}, E = 1.0, ν = 0.3) where {dim, T, iT, fT, CellType}
    # _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    # if _T <: Integer
    #     T = Float64
    # else
    #     T = _T
    # end
    if CellType === :Linear
        truss_grid = TrussGrid(node_points, elements, supports)
    else
        @assert false "not implemented"
    end

    # * load nodeset
    if haskey(truss_grid.grid.nodesets, "load")
        pop!(truss_grid.grid.nodesets, "load")
    end
    load_nodesets = Set{Int}()
    for (k,_) in loads
        push!(load_nodesets, k)
    end
    addnodeset!(truss_grid.grid, "load", load_nodesets)

    # * support nodeset
    for i=1:dim
        if haskey(truss_grid.grid.nodesets, "fixed_u$i")
            pop!(truss_grid.grid.nodesets, "fixed_u$i")
        end
        support_nodesets = Set{Int}()
        for kval in supports
            if kval[2][dim]
                push!(support_nodesets, kval[1])
            end
        end
        addnodeset!(truss_grid.grid, "fixed_u$i", support_nodesets)
    end

    # Create displacement field u
    dh = DofHandler(truss_grid.grid)
    if CellType === :Linear
        # truss linear
        ip = Lagrange{dim,RefCube,1}()
        # push!(dh, :u, dim, ip)
        push!(dh.field_names, :u)
        push!(dh.field_dims, dim)
        push!(dh.field_interpolations, ip)
        # push!(dh.bc_values, BCValues(ip, ip))
    else
        # TODO truss 2-order
        @assert false "not implemented"
        # ip = Lagrange{2, RefCube, 2}()
        # push!(dh, :u, dim, ip)
    end
    close!(dh)

    ch = ConstraintHandler(dh)
    #dbc1 = Dirichlet(:u, getfaceset(truss_grid.grid, "fixed_u1"), (x,t)->T[0], [1])
    # dbc2 = Dirichlet(:u, getnodeset(truss_grid.grid, "fixed_u2"), (x,t)->T[0], [2])
    for i=1:dim
        dbc = Dirichlet(:u, getnodeset(truss_grid.grid, "fixed_u$i"), (x,t)->T[0], [i])
        add!(ch, dbc)
    end
    close!(ch)

    # update the DBC to current time
    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)

    fnode = Tuple(getnodeset(truss_grid.grid, "load"))[1]
    node_dofs = metadata.node_dofs
    force_dof = node_dofs[2, fnode]

    # N = nnodespercell(truss_grid)
    # M = nfacespercell(truss_grid)

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)

    return TrussProblem(truss_grid, E, ν, ch, black, white, varind, metadata)
end

# nnodespercell(p::Union{PointLoadCantilever, HalfMBB}) = nnodespercell(p.truss_grid)
# function getcloaddict(p::Union{PointLoadCantilever{dim, T}, HalfMBB{dim, T}}) where {dim, T}
#     f = T[0, -p.force, 0]
#     fnode = Tuple(getnodeset(p.truss_grid.grid, "down_force"))[1]
#     return Dict{Int, Vector{T}}(fnode => f)
# end