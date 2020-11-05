using AbstractPlotting
using TopOpt.TopOptProblems: getdim
# TODO: build a MakieLayout

function draw_truss_problem!(scene, problem::TrussProblem;
    area=nothing, stress=nothing, line_width::Float64=1.0, 
    draw_supp::Bool=true, supp_scale::Float64=0.1,
    xaxis_label::String="x", plot_limits=undef)

    ndim = getdim(problem)
    ncells = JuAFEM.getncells(problem)

    if area !== nothing
        @assert(ncells == length(area))
        # a = reshape([area area]', 2*ncells)
        # a ./= maximum(a)
        # a .*= line_width
    else
        a = ones(ncells) * line_width
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

    # seg_ids = reshape([T[:,1] T[:,2]], length(T[:,1])+length(T[:,2]))
    # seg_ids = reshape(T', prod(size(T)))
    # if plot_limits == undef
    #     max_xlim = maximum(X[:,1]) - minimum(X[:,1])
    #     max_ylim = maximum(X[:,2]) - minimum(X[:,2])
    #     max_lim = max(max_xlim, max_ylim)
    #     plot_limits = FRect(minimum(X[:,1]), minimum(X[:,2]),
    #                    max_lim, max_lim)
    # end

    nodes = problem.truss_grid.grid.nodes
    PtT = ndim == 2 ? Point2f0 : Point3f0
    edges_pts = [PtT(nodes[cell.nodes[1]].x) => PtT(nodes[cell.nodes[2]].x) for cell in problem.truss_grid.grid.cells]
    @show edges_pts

    linesegments!(scene, edges_pts, 
                #   linewidth = a,
                  color = color,
                #   limits = plot_limits,
                  axis = (names = (axisnames = (xaxis_label, "y"),),
                          grid = (linewidth = (1, 1),),
                          )
                  )

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
