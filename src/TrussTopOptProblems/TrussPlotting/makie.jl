using AbstractPlotting: linesegments!, Point2f0, Point3f0
using AbstractPlotting.MakieLayout
using Makie
using TopOpt.TopOptProblems: getdim
using TrussTopOpt.TrussTopOptProblems: TrussProblem, get_fixities_node_set_name
using JuAFEM
using LinearAlgebra: norm

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
    else
        a = ones(2*ncells)
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

    if ndim == 2
        ax1 = layout[1, 1] = LAxis(scene)
        # tightlimits!(ax1)
        # ax1.aspect = AxisAspect(1)
        ax1.aspect = DataAspect()
    else
        # https://jkrumbiegel.github.io/MakieLayout.jl/v0.3/layoutables/#LScene-1
        # https://makie.juliaplots.org/stable/cameras.html#D-Camera
        ax1 = layout[1, 1] = LScene(scene, camera = cam3d!, raw = false)
    end
    # TODO show the ground mesh in another LAxis https://makie.juliaplots.org/stable/makielayout/grids.html
    # ax1.title = "Truss TopOpt result"

    # sl1 = layout[2, 1] = LSlider(scene, range = 0.01:0.01:10, startvalue = 1.0)
    lsgrid = labelslidergrid!(scene,
        ["support scale", "load scale", "arrow size", "vector linewidth", "element linewidth"],
        Ref(LinRange(0.01:0.01:10)); # same range for every slider via broadcast
        # formats = [x -> "$(round(x, digits = 2))$s" for s in ["", "", ""]],
        # width = 350,
        tellheight = false,
    )
    set_close_to!(lsgrid.sliders[1], 1.0)
    set_close_to!(lsgrid.sliders[2], 1.0)
    set_close_to!(lsgrid.sliders[3], 0.2)
    set_close_to!(lsgrid.sliders[4], 6.0)
    set_close_to!(lsgrid.sliders[5], 6.0)
    arrow_size = lift(s -> s, lsgrid.sliders[3].value)
    arrow_linewidth = lift(s -> s, lsgrid.sliders[4].value)
    layout[2, 1] = lsgrid.layout

    # * draw element
    # http://juliaplots.org/MakieReferenceImages/gallery//tutorial_linesegments/index.html
    element_linewidth = lift(s -> a.*s, lsgrid.sliders[5].value)
    linesegments!(ax1, edges_pts, 
                  linewidth = element_linewidth,
                  color = color)

    # fixties vectors
    for i=1:ndim
        nodeset_name = get_fixities_node_set_name(i)
        fixed_node_ids = JuAFEM.getnodeset(problem.truss_grid.grid, nodeset_name)
        dir = zeros(ndim)
        dir[i] = 1.0
        scaled_base_pts = lift(s->[PtT(nodes[node_id].x) - PtT(dir*s) for node_id in fixed_node_ids], 
            lsgrid.sliders[1].value)
        scaled_fix_dirs = lift(s->fill(PtT(dir*s), length(fixed_node_ids)), lsgrid.sliders[1].value)
        Makie.arrows!(
            ax1,
            scaled_base_pts,
            scaled_fix_dirs,
            arrowcolor=:orange,
            arrowsize=arrow_size,
            linecolor=:orange,
            linewidth=arrow_linewidth,
        )
    end
    # load vectors
    scaled_load_dirs = lift(s->[PtT(force/norm(force)*s) for force in values(problem.force)], 
        lsgrid.sliders[2].value)
    Makie.arrows!(
        ax1,
        [PtT(nodes[node_id].x) for node_id in keys(problem.force)],
        scaled_load_dirs,
        arrowcolor=:purple,
        arrowsize=arrow_size,
        linecolor=:purple,
        linewidth=arrow_linewidth,
    )
end
