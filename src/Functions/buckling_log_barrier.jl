# @params 
mutable struct BucklingLogBarrier{T} <: AbstractFunction{T}
	problem::StiffnessTopOptProblem
    solver::AbstractDisplacementSolver
    cheqfilter::AbstractCheqFilter
    comp::T # function value
    cell_comp::AbstractVector # value at each cell
    grad::AbstractVector
    tracing::Bool
    topopt_trace::TopOptTrace{T}
    reuse::Bool
    fevals::Integer
    maxfevals::Int
end
Utilities.getpenalty(c::BucklingLogBarrier) = c |> getsolver |> getpenalty
Utilities.setpenalty!(c::BucklingLogBarrier, p) = setpenalty!(getsolver(c), p)
TopOpt.dim(::BucklingLogBarrier) = 1

function BucklingLogBarrier(problem, solver::AbstractDisplacementSolver, args...; kwargs...)
    BucklingLogBarrier(problem, solver, args...; kwargs...)
end
function BucklingLogBarrier(problem::StiffnessTopOptProblem{dim, T}, solver::AbstractDisplacementSolver, ::Type{TI}=Int; rmin = T(0), filterT = nothing, tracing = false, logarithm = false, maxfevals = 10^8, preproj=nothing, postproj=nothing) where {dim, T, TI}

    rmin == 0 && filterT !== nothing && throw("Cannot use a filter radius of 0 in a density filter.")
    if filterT isa Nothing
        cheqfilter = SensFilter(Val(false), solver, rmin)
    elseif filterT isa SensFilter
        cheqfilter = SensFilter(Val(true), solver, rmin)
    else
        if preproj isa Nothing && postproj isa Nothing
            cheqfilter = DensityFilter(Val(true), solver, rmin)
        else
            cheqfilter = ProjectedDensityFilter(DensityFilter(Val(true), solver, rmin), preproj, postproj)
        end
    end
    comp = T(0)
    cell_comp = zeros(T, getncells(problem.ch.dh.grid))
    grad = fill(T(NaN), length(cell_comp) - sum(problem.black) - sum(problem.white))
    topopt_trace = TopOptTrace{T,TI}()
    reuse = false
    fevals = TI(0)
    return BucklingLogBarrier(problem, solver, cheqfilter, comp, cell_comp, grad, tracing, topopt_trace, reuse, fevals, maxfevals)
end

function (o::BucklingLogBarrier{T})(x, grad = o.grad) where {T}
    @unpack cell_comp, solver, tracing, cheqfilter, topopt_trace = o
    @unpack elementinfo, u, xmin = solver
    @unpack metadata, Kes, black, white, varind = elementinfo
    @unpack cell_dofs = metadata

    @timeit to "Eval obj and grad" begin
        penalty = getpenalty(o)
        if cheqfilter isa AbstractDensityFilter
            copyto!(solver.vars, cheqfilter(x))
        else
            copyto!(solver.vars, x)
        end
        if o.reuse
            if !tracing
                o.reuse = false
            end
        else
            o.fevals += 1
            setpenalty!(solver, penalty.p)
            solver()
        end
        obj = compute_compliance(cell_comp, grad, cell_dofs, Kes, u, 
                                 black, white, varind, solver.vars, penalty, xmin)

        o.comp = obj

        if cheqfilter isa AbstractDensityFilter
            grad .= TopOpt.jtvp!(similar(grad), cheqfilter, x, grad)
        else
            cheqfilter(grad)
        end
        if o.grad !== grad
            copyto!(o.grad, grad)
        end
        
        if o.tracing
            if o.reuse
                o.reuse = false
            else
                push!(topopt_trace.c_hist, obj)
                push!(topopt_trace.x_hist, copy(x))
                if length(topopt_trace.x_hist) == 1
                    push!(topopt_trace.add_hist, 0)
                    push!(topopt_trace.rem_hist, 0)
                else
                    push!(topopt_trace.add_hist, sum(topopt_trace.x_hist[end] .> topopt_trace.x_hist[end-1]))
                    push!(topopt_trace.rem_hist, sum(topopt_trace.x_hist[end] .< topopt_trace.x_hist[end-1]))
                end
            end
        end
    end
    return o.comp::T
end

function compute_buckling_log_barrier!(cell_comp::Vector{T}, grad, cell_dofs, Kes, u, 
                            black, white, varind, x, penalty, xmin) where {T}
    obj = zero(T)
    @inbounds for i in 1:size(cell_dofs, 2)
        cell_comp[i] = zero(T)
        Ke = rawmatrix(Kes[i])
        for w in 1:size(Ke,2)
            for v in 1:size(Ke, 1)
                cell_comp[i] += u[cell_dofs[v,i]]*Ke[v,w]*u[cell_dofs[w,i]]
            end
        end

        if black[i]
            obj += cell_comp[i]
        elseif white[i]
            if PENALTY_BEFORE_INTERPOLATION
                obj += xmin * cell_comp[i] 
            else
                p = penalty(xmin) * cell_comp[i]
            end
        else
            d = ForwardDiff.Dual{T}(x[varind[i]], one(T))
            if PENALTY_BEFORE_INTERPOLATION
                p = density(penalty(d), xmin)
            else
                p = penalty(density(d, xmin))
            end
            grad[varind[i]] = -p.partials[1] * cell_comp[i]
            obj += p.value * cell_comp[i]
        end
    end

    return obj    
end


