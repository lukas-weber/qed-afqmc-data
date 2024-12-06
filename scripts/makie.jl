using CairoMakie
using MakieCore

const marker_cycle =
    [:circle, :diamond, :rect, :utriangle, :ltriangle, :dtriangle, :rtriangle]

function qed_theme()
    return merge(
        theme_latexfonts(),
        Theme(
            figure_padding = 1,
            Axis = (; xgridvisible = false, ygridvisible = false),
            Legend = (;
                framevisible = false,
                labelsize = 12,
                titlesize = 14,
                patchsize = (14, 14),
            ),
            palette = (color = Makie.wong_colors(), marker = marker_cycle),
            MeasPlot = (cycle = Cycle([:color, :marker], covary = true),),
            Scatter = (cycle = Cycle([:color, :marker], covary = true),),
        ),
    )
end

@recipe MeasPlot (x, y) begin
    errorcolor = @inherit markercolor
    errorsinfront = false
    MakieCore.documented_attributes(Scatter)...
end

function Makie.plot!(mp::MeasPlot{<:Tuple{AbstractVector,AbstractVector}})
    vals = @lift(Measurements.value.($(mp.y)))
    errs = @lift(Measurements.uncertainty.($(mp.y)))
    if !mp.errorsinfront[]
        errorbars!(mp, mp.x, vals, errs, color = mp.errorcolor)
    end
    scatter!(
        mp,
        mp.x,
        vals;
        marker = mp.marker,
        color = mp.color,
        strokecolor = mp.strokecolor,
        strokewidth = mp.strokewidth,
        markersize = mp.markersize,
    )
    if mp.errorsinfront[]
        errorbars!(mp, mp.x, vals, errs, color = mp.errorcolor)
    end
    return mp
end

function Makie.legendelements(plot::MeasPlot, legend)
    les = LegendElement[
        LineElement(
            linepoints = [Point2f(0.5, 0.1), Point2f(0.5, 0.9)],
            color = plot.errorcolor,
        ),
        # MarkerElement(points = [Point2f(0.5, 0), Point2f(0.5, 1)], marker = :hline, markersize = 10),
        MarkerElement(
            points = [Point2f(0.5, 0.5)],
            marker = plot.marker,
            color = plot.color,
            strokecolor = plot.strokecolor,
            strokewidth = plot.strokewidth,
            markersize = plot.markersize,
        ),
    ]

    if plot.errorsinfront[]
        return reverse(les)
    end
    return les
end
