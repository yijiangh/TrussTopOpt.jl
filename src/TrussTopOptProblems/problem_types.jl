using TopOpt.TopOptProblems: StiffnessTopOptProblem, Metadata
# using TrussTopOpt.TrussTopOptProblems: getcelldim

# @params struct TrussProblem{xdim, T, N, M} <: StiffnessTopOptProblem{xdim, T}
struct TrussProblem{xdim,T,N,M} <: StiffnessTopOptProblem{xdim,T}
    truss_grid::TrussGrid{xdim,T,N,M}
    E::T
    ν::T
    ch::ConstraintHandler{<:DofHandler{xdim,<:JuAFEM.Cell{xdim,N,M},T},T}
    force::Dict{Int, SVector{xdim, T}}
    black::AbstractVector
    white::AbstractVector
    varind::AbstractVector{Int} # full dof -> free dof, based on black & white
    metadata::Metadata
end

function TrussProblem(::Type{Val{CellType}}, node_points::Dict{iT, SVector{xdim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
    loads::Dict{iT, SVector{xdim, T}}, supports::Dict{iT, SVector{xdim, fT}}, E = 1.0, ν = 0.3) where {xdim, T, iT, fT, CellType}
    # _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    # if _T <: Integer
    #     T = Float64
    # else
    #     T = _T
    # end
    if CellType === :Linear
        truss_grid = TrussGrid(node_points, elements, supports)
        # we assume geom interpolation order = function interpolation order
        geom_order = 1
    else
        @assert false "Other cell type not implemented"
    end
    # reference domain dimension
    ξdim = 1

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
    for i=1:xdim
        if haskey(truss_grid.grid.nodesets, "fixed_u$i")
            pop!(truss_grid.grid.nodesets, "fixed_u$i")
        end
        support_nodesets = Set{Int}()
        for kval in supports
            if kval[2][xdim]
                push!(support_nodesets, kval[1])
            end
        end
        addnodeset!(truss_grid.grid, "fixed_u$i", support_nodesets)
    end

    # Create displacement field u
    dh = DofHandler(truss_grid.grid)
    if CellType === :Linear
        # truss linear
        # interpolation_space
        ip = Lagrange{ξdim, RefCube, geom_order}()
        push!(dh, :u, xdim, ip)
    else
        # TODO truss 2-order
        @assert false "not implemented"
        # ip = Lagrange{2, RefCube, 2}()
        # push!(dh, :u, xdim, ip)
    end
    close!(dh)

    ch = ConstraintHandler(dh)
    # @show getnodeset(truss_grid.grid, "fixed_u1")
    # @show getnodeset(truss_grid.grid, "fixed_u2")
    for i=1:xdim
        dbc = Dirichlet(:u, getnodeset(truss_grid.grid, "fixed_u$i"), (x,t)->T[0], [i])

        # @show field_idx = JuAFEM.find_field(ch.dh, dbc.field_name)
        # @show interpolation = JuAFEM.getfieldinterpolation(ch.dh, field_idx)#ch.dh.field_interpolations[field_idx]
        # @show field_dim = JuAFEM.getfielddim(ch.dh, field_idx)#ch.dh.field_dims[field_idx]
        # @show bcvalue = JuAFEM.getbcvalue(ch.dh, field_idx)

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

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)

    return TrussProblem(truss_grid, E, ν, ch, loads, black, white, varind, metadata)
end

TopOpt.TopOptProblems.nnodespercell(p::TrussProblem) = nnodespercell(p.truss_grid)
# function getcloaddict(p::Union{PointLoadCantilever{dim, T}, HalfMBB{dim, T}}) where {dim, T}
#     f = T[0, -p.force, 0]
#     fnode = Tuple(getnodeset(p.truss_grid.grid, "down_force"))[1]
#     return Dict{Int, Vector{T}}(fnode => f)
# end

# https://github.com/lijas/JuAFEM.jl/blob/line2/src/Grid/grid.jl
const Line2d = Cell{2,2,2}
const Line3d = Cell{3,2,2}
const QuadraticLine = Cell{1,3,2}

# 1D: vertices
JuAFEM.faces(c::Union{Line,QuadraticLine}) = (c.nodes[1], c.nodes[2])
JuAFEM.vertices(c::Union{Line,Line2d,Line3d,QuadraticLine}) = (c.nodes[1], c.nodes[2])

# 2D: vertices, faces
JuAFEM.faces(c::Line2d) = ((c.nodes[1],c.nodes[2]),) 

# 3D: vertices, edges, faces
JuAFEM.edges(c::Line3d) = ((c.nodes[1],c.nodes[2]),) 

JuAFEM.default_interpolation(::Union{Type{Line},Type{Line2d},Type{Line3d}}) = Lagrange{1,RefCube,1}()
JuAFEM.default_interpolation(::Type{QuadraticLine}) = Lagrange{1,RefCube,2}()