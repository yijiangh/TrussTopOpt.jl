using TopOpt.TopOptProblems: AbstractGrid

const Vec = JuAFEM.Vec

# @params struct TrussGrid{xdim,N,M,C<:JuAFEM.Cell{xdim,N,M},T} <: AbstractGrid{xdim, T}
struct TrussGrid{xdim,T,N,M,TG<:JuAFEM.Grid{xdim,<:JuAFEM.Cell{xdim,N,M},T}} <: AbstractGrid{xdim, T}
    grid::TG
    white_cells::BitVector
    black_cells::BitVector
    constant_cells::BitVector
    crosssecs::Vector{T}
end
# nels::NTuple{dim, Int} # num of elements in x,y,z direction in the ground mesh, not applicable to truss
# sizes::NTuple{dim, T}  # length of the ground mesh in x,y,z direction, not applicaiton to truss
# corners::NTuple{2, Vec{dim, T}} # corner of the ground mesh

nnodespercell(::TrussGrid{xdim,T,N,M}) where {xdim,T,N,M} = N
nfacespercell(::TrussGrid{xdim,T,N,M}) where {xdim,T,N,M} = M
nnodes(cell::Type{JuAFEM.Cell{dim,N,M}}) where {dim, N, M} = N
nnodes(cell::JuAFEM.Cell) = nnodes(typeof(cell))

function TrussGrid(node_points::Dict{iT, SVector{xdim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
        boundary::Dict{iT, SVector{xdim, fT}}; crosssecs=undef) where {xdim, T, iT, fT}
    # ::Type{Val{CellType}}, 
    # if CellType === :Linear
    #     geoshape = Line
    # else
    #     @assert false "not implemented"
    #     # geoshape = QuadraticQuadrilateral
    # end

    grid = _LinearTrussGrid(node_points, elements, boundary)
    ncells = getncells(grid)
    if crosssecs == undef
        # default cross section area to 1
        crosssecs = ones(T, ncells)
    end
    return TrussGrid(grid, falses(ncells), falses(ncells), falses(ncells), crosssecs)
end

function _LinearTrussGrid(node_points::Dict{iT, SVector{xdim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
        boundary::Dict{iT, SVector{xdim, fT}}) where {xdim, T, iT, fT}
    n_nodes = length(node_points)

    # * Generate cells, Line2d or Line3d
    CellType = Cell{xdim,2,2}
    cells = CellType[]
    for e in elements
        push!(cells, CellType((e[2]...,)))
    end

    # * Generate nodes
    nodes = Node{xdim,T}[]
    for kval in node_points
        push!(nodes, Node((kval[2]...,)))
    end

    # * label boundary (cell, face)
    cell_from_node = node_neighbors(cells)
    boundary_conditions = Tuple{Int,Int}[]
    for (v, _) in boundary
        for c in cell_from_node[v]
            push!(boundary_conditions, (c, v))
        end
    end
    boundary_matrix = JuAFEM.boundaries_to_sparse(boundary_conditions)

    # # Cell face sets
    # facesets = Dict("left"  => Set{Tuple{Int,Int}}([boundary[1]]),
    #                 "right" => Set{Tuple{Int,Int}}([boundary[2]]))

    return Grid(cells, nodes, boundary_matrix=boundary_matrix)
end

"""
Compute a dict for querying connected cells (elements) to a given node index

Returns: dict (node_idx => [cell_idx, ...])
"""
function node_neighbors(cells)
    # Contains the elements that each node contain
    cell_from_node = Dict{Int, Set{Int}}()
    for (cellid, cell) in enumerate(cells)
        for v in cell.nodes
            if !haskey(cell_from_node, v)
                cell_from_node[v] = Set{Int}()
            end
            push!(cell_from_node[v], cellid)
        end
    end
    cell_from_node
end