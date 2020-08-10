using TopOpt.TopOptProblems: AbstractGrid

const Vec = JuAFEM.Vec

struct TrussGrid{dim, T, N, M, TG<:JuAFEM.Grid{dim, N, T, M}} <: AbstractGrid{dim, T}
    grid::TG
    white_cells::BitVector
    black_cells::BitVector
    constant_cells::BitVector
end
    # nels::NTuple{dim, Int}
    # sizes::NTuple{dim, T}
    # corners::NTuple{2, Vec{dim, T}}

nnodespercell(::TrussGrid{dim,T,N,M}) where {dim, T, N, M} = N
nfacespercell(::TrussGrid{dim,T,N,M}) where {dim, T, N, M} = M

nnodes(cell::Type{JuAFEM.Cell{dim,N,M}}) where {dim, N, M} = N
nnodes(cell::JuAFEM.Cell) = nnodes(typeof(cell))

function TrussGrid(node_points::Dict{iT, SVector{dim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
        boundary::Dict{iT, SVector{dim, fT}}) where {dim, T, iT, fT}
    # ::Type{Val{CellType}}, 
    # if CellType === :Linear
    #     geoshape = Line
    # else
    #     @assert false "not implemented"
    #     # geoshape = QuadraticQuadrilateral
    # end

    grid = _LinearTrussGrid(node_points, elements, boundary)
    ncells = getncells(grid)
    return TrussGrid(grid, falses(ncells), falses(ncells), falses(ncells))
end

function _LinearTrussGrid(node_points::Dict{iT, SVector{dim, T}}, elements::Dict{iT, Tuple{iT, iT}}, 
        boundary::Dict{iT, SVector{dim, fT}}) where {dim, T, iT, fT}
    n_nodes = length(node_points)
    TrussLine = Cell{dim,2,2}
    # TrussLine = Line

    # * Generate cells
    cells = TrussLine[]
    for e in elements
        push!(cells, TrussLine((e[2]...,)))
    end

    # * Generate nodes
    nodes = Node{dim,T}[]
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