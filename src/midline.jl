using LinearAlgebra
using Statistics

function midline_points(img, kernel = KernelFactors.ando5, q=0.995)
    ry, rx = real.(imgradients(img, kernel))
    ryy, ryx = imgradients(ry, kernel)
    rxy, rxx = imgradients(rx, kernel)
    hessians = broadcast( (xx,xy,yy)->Symmetric([xx xy; xy yy]), rxx, rxy, ryy )
    e = eigen.(hessians)
    imax = [argmax(abs.(x.values)) for x in e]
    val_max = [x.values[i] for (x,i) in zip(e,imax)]
    vec_max = [x.vectors[:,i] for (x,i) in zip(e,imax)]
    nx = getindex.(vec_max,1)
    ny = getindex.(vec_max,2)
    t = @. - (rx*nx + ry*ny) / (rxx*nx^2 + 2rxy*nx*ny + ryy*ny^2)
    px = t .* nx
    py = t .* ny
    th = quantile(filter(!isnan,abs.(vec(val_max))), q)
    #th = mean(extrema(abs.(val2)))

    i = @. (-0.5 < px < 0.5) & (-0.5 < py < 0.5) & (abs(val_max) > th)
    centers = [(i[2] + px[i], i[1] + py[i]) for i in findall(i)]
    centers, i, px, py, val_max
end

function link_points(mask, val2, nx, ny)
    i = findall(mask)
    i_max = i[argmax(val2[i])]
    result = Set{CartesianIndex{2}}((i_max,))
    grow_line!(result, i_max, nx, ny)

end

function grow_line!(result, i, nx, ny)
    # two candidate directions
    d1 = [-ny[i], nx[i]]
    d2 = [ny[i], -nx[i]]
    #...
end


##

using ImageFiltering
using GeometryTypes
centered = ImageFiltering.centered

function line_ends(thin_mask)
    kern = kernelfactors(centered.(([1,1,1],[1,1,1])))
    findall(thin_mask .& (imfilter(thin_mask, kern) .== 2))
end

const _neighborhood_car = CartesianIndex.([(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1,),(0,-1)])
const _neighborhood = Point2.([(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1,),(0,-1)])

_ind2point2(x::CartesianIndex) = Point2(x[2],x[1])
function follow_line(thin_mask)
    mask_indices = CartesianIndices(thin_mask)
    ends = line_ends(thin_mask)
    isempty(ends) && error("No ends found")

    done = false
    visited = Set{CartesianIndex{2}}()
    line = Point2{Int}[]
    curr = ends[1]
    while !done
        push!(line, _ind2point2(curr))
        push!(visited, curr)
        done = true
        for d in _neighborhood_car
            next = curr + d
            next in mask_indices || continue
            if thin_mask[next] && !(next in visited)
                done = false
                curr = next
                break
            end
        end
    end
    # if length(ends) < 2
    #     @warn("One-ended midline")
    # elseif curr != ends[2]
    #     @warn("Other end not reached")
    # end
    line
end

using ImageMorphology

default_threshold(image) = ((min,max) = extrema(image); min + 0.5*(max-min))

function midline_thinmask(image, th)
    mask = image .< th
    thinned = thinning(mask)
end


midline_by_thinning(image, th) = follow_line(midline_thinmask(image,th))

anglevec_by_thinning(image, s=0:0.025:1, g=0,
                     th=default_threshold(image)) = spline_anglevec(
                line2spline(imfilter(midline_by_thinning(image,th),Kernel.gaussian((g,)))), s)

import DSP

function spline_anglevec(spl, s=0:0.025:1)
    x, y = collect(eachrow(spl(s)))
    dx, dy = diff(x), diff(y)
    θ = DSP.unwrap(atan.(dy,dx))
#    θ .-= mean(θ)
end

## Line caching

# const VideoID = NamedTuple{(:path, :idx),Tuple{String,Int64}}
#
# function frame_cache()
#     cache = Dict()
#     (path, idx) -> get!(cache, (path=path, idx=idx)) do
#         Elegans.read_video(Elegans.videopath_f(path),idx)
#     end
# end

summarize_exception( e, trace ) = sprint(e, trace) do io, e, trace
    showerror(io, e)
    println(io)
    show(io, MIME("text/plain"), stacktrace(trace))
end

# exc_f defaults to typeof, to keep only the type of exception:
#  - BoundError-s on images keep a reference to the image,
#    which makes them very large
#  - Converting such BoundError-s directly to string is extremely
#    slow (~1min!), so not a good option either
function trying_cache(f, ::Type{IN}, ::Type{OUT};
                      cache = Dict{IN,Union{OUT,DataType}}(),
                      exc_f = (e,t)->typeof(e)) where {IN, OUT}
    x -> get!(cache, x) do
        try
            f(x)
        catch e
            exc_f( e, catch_backtrace() )
        end
    end
end

using GeometryTypes


function midline_cache(videocache, threshold, kernel = Kernel.gaussian(1))
    trying_cache( Int, Vector{Point2{Int}} ) do idx
        frame = get_frame(videocache,idx)
        midline_by_thinning(imfilter(Gray.(frame)[1:end-1,1:end-1], kernel), threshold)
    end
end

using Dierckx

spline_cache(midlines) = trying_cache( i->line2spline(midlines(i)),
                                       Int, ParametricSpline )

# function midline_cache(videocache, kernel = Kernel.gaussian(1))
#     #cache = Dict{ Int, Union{ Vector{Point{2,Int}}, Exception } }()
#     cache = Dict{ Int, Union{ Vector{Point{2,Int}}, DataType } }()
#     idx -> get!(cache, idx) do
#         try
#             frame = get_frame(videocache,idx)
#             midline_by_thinning(imfilter(Gray.(frame)[1:end-1,1:end-1], kernel))
#         catch e
#             typeof(e)
#         end
#     end
# end
#
# function spline_cache(midlines)
#     cache = Dict{ Int, Union{ ParametricSpline, DataType } }()
#     idx -> get!(cache, idx) do
#         try
#             line2spline(midlines(idx))
#         catch e
#             typeof(e)
#         end
#     end
# end
