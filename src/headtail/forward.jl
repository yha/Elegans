using ImageFiltering.KernelFactors: gaussian
using GeometryBasics, LinearAlgebra
using Statistics
using StatsBase
using ProgressLogging: @progress

function forward_speed(x, y, d, k)
    @assert length(x) == length(y) == length(d) + 1

    x = replace(x, missing => NaN)
    y = replace(y, missing => NaN)
    dx = diff(imfilter(x, k, NA()))
    dy = diff(imfilter(y, k, NA()))

    i = .!ismissing.(d)

    v = Point2.(dx, dy)
    vd = dot.(v[i], normalize.(d[i]))
end

kern(σ) = gaussian(σ, 6 * ceil(Int, σ) + 1)
speed_heatmap_data(σ, x, y, d) =
    @progress fwd = [filter(!isnan, forward_speed(x, y, d, kern(σ))) for σ in σ]
speed_heatmap_plot((from,to)::Tuple{Integer,Integer}, σ, fwd) = _speed_heatmap(range(from, to, length=100), σ, fwd)
function speed_heatmap_plot(speed_edges, σ, fwd)
    speeds = midpoints(speed_edges)
    heatmap(speeds, σ,
        reduce(hcat, fit(Histogram, f, speed_edges).weights for f in fwd)',
        #c=:grays_r,
        xlabel = "speed",
        #ylabel = "\$\\sigma\$",
        ylabel = "time",
        legend = false,
        colorbar = true,
    )
    vline!([0], c = "white", ls = :dot)
    plot!(mean.(fwd), σ)
end
function speed_heatmap(speed_edges, σ, x, y, d)
    fwd = speed_heatmap_data(σ, x, y, d)
    speed_heatmap_plot(speed_edges,σ,fwd)
end
