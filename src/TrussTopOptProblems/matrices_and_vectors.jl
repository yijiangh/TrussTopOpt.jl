using TopOpt.TopOptProblems: getgeomorder, getdh, getE, getν, getdensity, gettypes
using LinearAlgebra: norm
import JuAFEM: getngeobasefunctions, getn_scalarbasefunctions

"""
Generate element stiffness matrices
"""
function make_Kes_and_fes(problem::TrussProblem, quad_order=1)
    make_Kes_and_fes(problem, quad_order, Val{:Static})
end

function make_Kes_and_fes(problem::TrussProblem, ::Type{Val{mat_type}}) where mat_type
    make_Kes_and_fes(problem, 1, Val{mat_type})
end

function make_Kes_and_fes(problem::TrussProblem{xdim, T}, quad_order, ::Type{Val{mat_type}}) where {xdim, T, mat_type}
    geom_order = getgeomorder(problem)
    dh = getdh(problem)
    E = getE(problem)
    # ν = getν(problem)
    # ρ = getdensity(problem)

    refshape = JuAFEM.getrefshape(dh.field_interpolations[1])

    # Shape functions and quadrature rule
    interpolation_space = Lagrange{ξdim, refshape, geom_order}()
    quadrature_rule = QuadratureRule{ξdim, refshape}(quad_order)
    cellvalues = CellScalarValues(T, quadrature_rule, interpolation_space; xdim=xdim)

    # Calculate element stiffness matrices
    n_basefuncs = getnbasefunctions(cellvalues)
    Kesize = xdim*n_basefuncs
    MatrixType, VectorType = gettypes(T, Val{mat_type}, Val{Kesize})
    Kes, weights = _make_Kes_and_weights(dh, Tuple{MatrixType, VectorType}, Val{n_basefuncs}, Val{xdim*n_basefuncs}, E, quadrature_rule, cellvalues)
    # dloads = _make_dloads(weights, problem, facevalues)

    # return Kes, weights, dloads, cellvalues, facevalues
    return Kes, weights, cellvalues
end

"""
    Kes, weights = _make_Kes_and_weights(dof_handler, Tuple{MatrixType, VectorType}, Val{n_basefuncs}, Val{dim*n_basefuncs}, E, quadrature_rule, cellvalues)
"""
function _make_Kes_and_weights(
    dh::DofHandler{xdim, N, T},
    ::Type{Tuple{MatrixType, VectorType}},
    ::Type{Val{n_basefuncs}},
    ::Type{Val{Kesize}},
    E::T, quadrature_rule, cellvalues::CellScalarValues{1,xdim}) where {xdim, N, T, MatrixType <: StaticArray, VectorType, n_basefuncs, Kesize}
    nel = getncells(dh.grid)
    Kes = Symmetric{T, MatrixType}[]
    sizehint!(Kes, nel)
    weights = [zeros(VectorType) for i in 1:nel]
    Ke_e = zeros(T, xdim, xdim)
    fe = zeros(T, Kesize)
    Ke_0 = Matrix{T}(undef, Kesize, Kesize)

    celliterator = CellIterator(dh)
    for (k, cell) in enumerate(celliterator)
        Ke_0 .= 0
        truss_reinit!(cellvalues, cell)
        fe = weights[k]
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for b in 1:n_basefuncs
                ∇ϕb = shape_gradient(cellvalues, q_point, b)
                ϕb = shape_value(cellvalues, q_point, b)
                for d2 in 1:xdim
                    # fe = @set fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
                    for a in 1:n_basefuncs
                        ∇ϕa = shape_gradient(cellvalues, q_point, a)
                        Ke_e .= E * ∇ϕa ⊗ ∇ϕb * dΩ
                        for d1 in 1:xdim
                            #if dim*(b-1) + d2 >= dim*(a-1) + d1
                            Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                            #end
                        end
                    end
                end
            end
        end
        weights[k] = fe
        if MatrixType <: SizedMatrix # Work around because full constructor errors
            push!(Kes, Symmetric(SizedMatrix{Kesize,Kesize,T}(Ke_0)))
        else
            push!(Kes, Symmetric(MatrixType(Ke_0)))
        end
    end
    return Kes, weights
end

function truss_reinit!(ci::CellIterator{xdim,C}, i::Int) where {xdim,C}
    ci.current_cellid[] = ci.cellset[i]

    ci.flags.nodes  && cellnodes!(ci.nodes, ci.dh, ci.current_cellid[])
    ci.flags.coords && cellcoords!(ci.coords, ci.dh, ci.current_cellid[])

    if isdefined(ci, :dh) && ci.flags.celldofs # update celldofs
        celldofs!(ci.celldofs, ci.dh, ci.current_cellid[])
    end
    return ci
end

@inline function truss_reinit!(cv::CellValues{ξdim,xdim,T}, ci::CellIterator{xdim,N,T}) where {ξdim,xdim,N,T}
    check_compatible_geointerpolation(cv, ci)
    reinit!(cv, ci.coords)
end

"""
Reinit a cell for a truss element, using the nodal coordinates `x`, cross section `A`
"""
function truss_reinit!(cv::CellValues{1,xdim}, x::AbstractVector{Vec{xdim,T}}, A::T) where {xdim,T}
    n_geom_basefuncs = getngeobasefunctions(cv)
    n_func_basefuncs = getn_scalarbasefunctions(cv)
    @assert length(x) == n_geom_basefuncs
    isa(cv, CellVectorValues) && (n_func_basefuncs *= xdim)

    @inbounds for i in 1:length(cv.qr_weights)
        w = cv.qr_weights[i]
        dxdξ = zero(Tensor{1,xdim})
        for j in 1:n_geom_basefuncs
            # in a truss element, x_j ∈ R, dMdξ_j ∈ R, ξ ∈ R
            # cv.dMdξ[j, i] is a 1-1 tensor here
            dxdξ += x[j] * cv.dMdξ[j, i][1]
        end
        # detJ = √(J' J), J = dxdξ ∈ R(n x 1)
        detJ = norm(dxdξ)
        detJ > 0.0 || throw(ArgumentError("det(J) is not positive: det(J) = $(detJ)"))
        cv.detJdV[i] = detJ * w * A
        Jinv = pinv(dxdξ)
        for j in 1:n_func_basefuncs
            # cv.dNdξ[j, i] is a 1-1 tensor here
            cv.dNdx[j, i] = cv.dNdξ[j, i][1] * Jinv'
        end
    end
end

