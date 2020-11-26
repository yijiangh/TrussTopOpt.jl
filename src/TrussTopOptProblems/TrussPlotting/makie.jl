using AbstractPlotting: linesegments!, Point2f0, Point3f0
using AbstractPlotting.MakieLayout
using TopOpt.TopOptProblems: getdim
using TrussTopOpt.TrussTopOptProblems: TrussProblem
using JuAFEM

"""
    scene, layout = draw_truss_problem(problem; crosssecs=result.topology)
"""
function draw_truss_problem(problem::TrussProblem; kwargs...)
    scene, layout = layoutscene() #resolution = (1200, 900)
    draw_truss_problem!(scene, layout, problem; kwargs...)
    display(scene)
    return scene, layout
end

function draw_truss_problem!(scene, layout, problem::TrussProblem;
    crosssecs=nothing, stress=nothing, linewidth::Float64=6.0)
    ndim = getdim(problem)
    ncells = JuAFEM.getncells(problem)

    if crosssecs !== nothing
        @assert(ncells == length(crosssecs))
        a = reshape([crosssecs crosssecs]', 2*ncells)
        # a ./= maximum(a)
        a .*= linewidth
    else
        a = ones(ncells) * linewidth
    end
    if stress !== nothing
        @assert(ncells == length(stress))
        q_color = Array{RGBAf0, 1}(undef, length(stress))
        for i=1:ncells
            if stress[i] < 0
                q_color[i] = RGBAf0(0,0,1,0.8)
            else
                q_color[i] = RGBAf0(1,0,0,0.8)
            end
        end
        color = reshape([q_color q_color]', 2*length(stress))
    else
        color = :black
    end

    nodes = problem.truss_grid.grid.nodes
    PtT = ndim == 2 ? Point2f0 : Point3f0
    edges_pts = [PtT(nodes[cell.nodes[1]].x) => PtT(nodes[cell.nodes[2]].x) for cell in problem.truss_grid.grid.cells]

    # TODO show the ground mesh in another LAxis https://makie.juliaplots.org/stable/makielayout/grids.html
    ax1 = layout[1, 1] = LAxis(scene) #, title = "Axis 1")
    tightlimits!(ax1)
    # ax1.title = "Truss TopOpt result"
    # ax1.aspect = AxisAspect(1)
    ax1.aspect = DataAspect()

    # http://juliaplots.org/MakieReferenceImages/gallery//tutorial_linesegments/index.html
    linesegments!(ax1, edges_pts, 
                  linewidth = a,
                  color = color)

    # if ndim == 2
    #     axis = scene[Axis2D]
    # else
    #     axis = scene[Axis3D]
    # end
    # axis[:scale] = (1.0, 1.0, 1.0)

    # if draw_supp
    #     fix_ids = findall(x->x==1, S)
    #     for s in fix_ids
    #         if 1 == s[2]
    #             # x dir
    #             arrows!(scene, [X[s[1], 1]], [X[s[1], 2]], [supp_scale], [0], linecolor = :green, arrowcolor = :green, limits = plot_limits,
    #             # linewidth = supp_scale*10,
    #             arrowsize = 0.1,
    #             )
    #         end
    #         if 2 == s[2]
    #             # y dir
    #             arrows!(scene, [X[s[1], 1]], [X[s[1], 2]], [0], [supp_scale], linecolor = :green, arrowcolor = :green, limits = plot_limits,
    #             # linewidth = supp_scale*10,
    #             arrowsize = 0.1,
    #             )
    #         end
    #         # if 1 == S[i,4]
    #         #     scatter!(scene, [X[S[i,1], 1]], [X[S[i,1], 2]], color = :green, markersize=supp_scale, marker=:x, limits = plot_limits)
    #         # else
    #         #     scatter!(scene, [X[S[i,1], 1]], [X[S[i,1], 2]], color = :green, markersize=supp_scale, marker=:circle, limits = plot_limits)
    #         # end
    #     end
    # end
    # axis = scene[Axis]
    # axis[:grid][:linewidth] = (1, 1)
end

function draw_load!(scene, X::Matrix{Float64}, F::Matrix{Float64}; load_scale::Float64=0.1, xaxis_label::String="x", plot_limits=undef)
    if plot_limits == undef
        max_xlim = maximum(X[:,1]) - minimum(X[:,1])
        max_ylim = maximum(X[:,2]) - minimum(X[:,2])
        max_lim = max(max_xlim, max_ylim)
        plot_limits = FRect(minimum(X[:,1]), minimum(X[:,2]),
                       max_lim, max_lim)
    end

    for i=1:size(F,1)
        if any(i -> i!=0, F[i,:])
            arrows!(scene, [X[i, 1]],
                           [X[i, 2]],
                           [F[i,1]*load_scale],
                           [F[i,2]*load_scale],
                           linecolor = :orange, arrowcolor = :orange,
                           # linewidth = load_scale*10,
                           arrowsize = 0.1,
                           limits = plot_limits,
                           axis = (names = (axisnames = (xaxis_label, "y"),),)
                           )
         end
    end
    # axis = scene[Axis]
    # axis[:grid][:linewidth] = (1, 1)
end
