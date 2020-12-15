# modified from https://github.com/mohamed82008/LinearElasticity.jl
using Einsum: @einsum
using LinearAlgebra: I

function get_Kσs(problem::TrussProblem{xdim, TT}, dofs, cellvalues) where {xdim, TT}
    Es = getE(problem)
    νs = getν(problem)
    As = getA(problem)

    dh = problem.ch.dh
    n = ndofs_per_cell(dh)
    global_dofs = zeros(Int, n)
    # element geometric stiffness matrix
    Kσs = [zeros(TT, n, n) for i in 1:getncells(dh.grid)]
    Kσ_e = zeros(TT, n, n)
    ψ_e = zeros(TT, xdim*n, xdim*n)
    G = zeros(TT, xdim*n, n)
    δ = Matrix{TT}(I, xdim, xdim)
    ϵ = zeros(TT, xdim, xdim)
    σ = zeros(TT, xdim, xdim)
    u = zeros(TT, xdim, xdim)

    n_basefuncs = getnbasefunctions(cellvalues)
    for (cellidx, cell) in enumerate(CellIterator(dh))
        Kσ_e .= 0
        truss_reinit!(cellvalues, cell, As[cellidx])
        celldofs!(global_dofs, dh, cellidx)
        E = Es[cellidx]
        ν = νs[cellidx]
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for d in 1:xdim
                ψ_e[(d-1)*xdim+1:d*xdim, (d-1)*xdim+1:d*xdim] .= 0
            end
            for a in 1:n_basefuncs
                ∇ϕ = shape_gradient(cellvalues, q_point, a)
                _u = @view dofs[(@view global_dofs[xdim*(a-1) .+ (1:xdim)])]
                @einsum u[i,j] = _u[i]*∇ϕ[j]
                @einsum ϵ[i,j] = 1/2*(u[i,j] + u[j,i])
                @einsum σ[i,j] = E*ν/(1-ν^2)*δ[i,j]*ϵ[k,k] + E*ν*(1+ν)*ϵ[i,j]
                for d in 1:xdim
                    ψ_e[(d-1)*xdim+1:d*xdim, (d-1)*xdim+1:d*xdim] .+= σ
                    G[(xdim*(d-1)+1):(xdim*d), (a-1)*xdim+d] .= ∇ϕ
                end
            end
            Kσ_e .+= G'*ψ_e*G*dΩ
        end
        Kσs[cellidx] .= Kσ_e
    end

    return Kσs
end

function buckling(problem::TrussProblem{xdim, T}, ginfo, einfo) where {xdim, T}
    dh = problem.ch.dh

    u = ginfo.K \ ginfo.f
    @show Kσs = get_Kσs(problem, u, einfo.cellvalues)
    Kσ = deepcopy(ginfo.K)

    if Kσ isa Symmetric
        Kσ.data.nzval .= 0
        assembler = JuAFEM.AssemblerSparsityPattern(Kσ.data, T[], Int[], Int[])
    else
        Kσ.nzval .= 0
        assembler = JuAFEM.AssemblerSparsityPattern(Kσ, T[], Int[], Int[])
    end

    global_dofs = zeros(Int, ndofs_per_cell(dh))
    Kσ_e = zeros(T, size(Kσs[1]))
    celliteratortype = CellIterator{typeof(dh).parameters...}
    _celliterator::celliteratortype = CellIterator(dh)
    TK = eltype(Kσs)
    for (i,cell) in enumerate(_celliterator)
        celldofs!(global_dofs, dh, i)
        if TK <: Symmetric
            JuAFEM.assemble!(assembler, global_dofs, Kσs[i].data)
        else
            JuAFEM.assemble!(assembler, global_dofs, Kσs[i])
        end
    end

    return ginfo.K, Kσ
end
