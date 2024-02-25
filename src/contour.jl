using Contour
import Images
using ImageFiltering
using Images: otsu_threshold, imfilter
using GeometryBasics: Point
using CircularArrays
using LinearAlgebra
using PolygonOps
using OffsetArrays: no_offset_view
import IterTools
import ShiftedArrays
const cshift = ShiftedArrays.circshift

# A closed 2d curve with circular indexing
struct Closed2DCurve{T} <: AbstractVector{Point{2,T}}
    vertices::CircularVector{Point{2,T}}
    function Closed2DCurve{T}(v::AbstractVector{Point{2,T}}; orient_pos=false) where T
        vertices = if orient_pos && !_is_positively_oriented(v)
            v = reverse(v)
            _is_positively_oriented(v) || @warn("Failed to ensure positive orientation. Curve might be ill-conditioned.")
            v
        else
            v
        end
        new{T}(CircularVector(vertices))
    end
end

Closed2DCurve(v::AbstractVector{Point{2,T}}; orient_pos=false) where T = Closed2DCurve{T}(v;orient_pos=orient_pos)

function Closed2DCurve( c::Contour.Curve2; orient_pos=false )
    _is_closed(c) || throw(ArgumentError("Curve not closed"))
    # Contour.Curve2 repeats the first vertex when closed.
    # Remove the repeated vertex.
    Closed2DCurve( Point.(c.vertices[1:end-1]); orient_pos=orient_pos )
end
vertices(c::Closed2DCurve) = c.vertices
vertices_closed(c::Closed2DCurve) = c.vertices[0:length(c.vertices)]

Base.size(c::Closed2DCurve) = (length(c.vertices),)
Base.getindex(c::Closed2DCurve, i::Int) = c.vertices[i]

Base.iterate(c::Closed2DCurve, state=firstindex(c.vertices)) = ( state > lastindex(c.vertices)
                ? nothing : (c.vertices[state], state+1) )

# This method is necessary to allow general indexing (e.g. with ranges) out of
# the underlying vector's axes.
Base.checkbounds(::Type{Bool}, ::Closed2DCurve, I...) = true

_circshift( c::Closed2DCurve, shift ) = circshift( c.vertices, shift )
Base.circshift( c::Closed2DCurve, shift ) = Closed2DCurve(_circshift(c))
function Base.circshift!( c::Closed2DCurve, shift )
    c.vertices .= _circshift( c, shift )
    c
end

_imfilter(c::Closed2DCurve, kernel) = imfilter( c.vertices, kernel, "circular" )
Images.imfilter(c::Closed2DCurve, kernel) = Closed2DCurve(_imfilter(c,kernel))


fwddiff(a) = a .- cshift(a,-1)
backdiff(a) = a .- cshift(a,1)
backmean(a) = (a .+ cshift(a,1)) ./ 2

function _curve_stats(vertices)
    dv = fwddiff(vertices)
    lens = norm.(dv)
    lens_mid = backmean(lens) # distances between edge midpoints
    θ = [atan(p[2],p[1]) for p in dv]
    dθ = fix_angle.(backdiff(θ))
    lens, lens_mid, dθ
end
function _total_curvature(vertices)
    _,_,dθ = _curve_stats(vertices)
    sum(dθ)
end
_turning_number(vertices) = _total_curvature(vertices) / 2π
_is_positively_oriented(vertices) = _total_curvature(vertices) > 0

function curvature( c::Closed2DCurve; warn_threshold=1e-3 )
    lens, lens_mid, dθ = _curve_stats(vertices(c))
    if any( lens_mid .< warn_threshold )
        @warn("Some curve points are very close ($(minimum(lens_mid)) pixels " *
              "between mid-segments along curve)", maxlog=4)
    end
    κ = dθ ./ lens_mid
end

_is_closed(c::Curve2) = c.vertices[end] == c.vertices[1]

function contour_lines( img, level )
    y, x = float.(axes(img))
    curves = lines(Contour.contour(x, y, img', level))
end

function closed_contour_levels(img::AbstractArray{<:AbstractFloat}, level::Number)::Vector{Closed2DCurve{Float64}}
    # Contour.contour does not support offset arrays, so keep and strip offsets
    # TODO: remove when offset arrays are directly supported
    offsets = first.(axes(img)) .-1
    img = no_offset_view(img)
    # Contour.contour(x,y,z) expects z to be indexed as
    # z[x[i],y[i]], contrary to the Images.jl convention img[y,x].
    # To get coordinates ordered as (x,y), we transpose img
    y, x = float.(axes(img))
    curves = lines(Contour.contour(x, y, img', level))
    curves = filter(_is_closed, curves)
    closed_curves = Closed2DCurve.(curves; orient_pos=true)
    closed_curves = [Closed2DCurve(c .+ Point2(offsets)) for c in closed_curves]
    sort( closed_curves; by=length, rev=true )
end
closed_contour_levels(img::AbstractArray{<:Number}, level::Number) = closed_contour_levels(float.(img), level)
closed_contour_levels(img::AbstractArray{<:Gray}, level) = closed_contour_levels(real.(img), real(level))

abstract type ContouringMethod end
struct Thresholding{T} <: ContouringMethod
    σ::Float64
    level::T
end
struct SeededSegmentation <: ContouringMethod
    σ::Float64
end

iscentered(a) = all(first(i) == -last(i) for i in axes(a))
assertcentered(a) = (@assert iscentered(a); a)

raw_worm_contours( img, method::Thresholding ) = closed_contour_levels(
                            imfilter( assertcentered(img), Kernel.gaussian(method.σ) ), method.level )
raw_worm_contours( img, method::SeededSegmentation ) = bgworm_segment_closed_contours(img, method.σ)

worm_contour( img, method ) = raw_worm_contours( img, method )[1]

incurve( p, c::Closed2DCurve; on=true ) = inpolygon( p, vertices_closed(c); in=true, on=on, out=false )

euclid_d(p,q) = norm(p-q)
# TODO optimize
incurve( c1::Closed2DCurve, c2::Closed2DCurve ) = all( incurve(p, c2) for p in c1 )
nearest_vertex( p, c::Closed2DCurve, d=euclid_d ) = argmin( Dict(q=>d(p,q) for q in c) )
nearest_vertex_index( p, c::Closed2DCurve, d=euclid_d ) = argmin( Dict(i=>d(p,q) for (i,q) in enumerate(c)) )
min_vertex_dist( p, c::Closed2DCurve, d=euclid_d ) = minimum( d(p,q) for q in c )


##
# Distance-based shape cutting (joining of inner and outer curve)

function cutlines(outer, inner, open, d=euclid_d)
    nearest_in_i = [nearest_vertex_index(p,inner,d) for p in open]
    nearest_out_i = [nearest_vertex_index(p,outer,d) for p in open]
    nearest_in = inner.vertices[nearest_in_i]
    nearest_out = outer.vertices[nearest_out_i]
    dist_in = [d(p,q) for (p,q) in zip(open,nearest_in)]
    dist_out = [d(p,q) for (p,q) in zip(open,nearest_out)]
    ins = dist_in .< dist_out
    cut_in = findfirst([x+1==y for (x,y) in IterTools.partition(ins,2,1)])
    cut_out = findfirst([x==y+1 for (x,y) in IterTools.partition(ins,2,1)])
    first_cut = cut_in === nothing ? nothing : (nearest_out_i[cut_in], nearest_in_i[cut_in+1], cut_in)
    second_cut = cut_out === nothing ? nothing : (nearest_out_i[cut_out+1], nearest_in_i[cut_out], cut_out)
    (first_cut, second_cut)
end

function cutline( outer, inner, open, d=euclid_d, σκ=2 )
    cut1, cut2 = cutlines( outer, inner, open, d )
    select_cutline( cut1, cut2, open, σκ )
end
function select_cutline( cutline1, cutline2, open, σκ=2 )
    κ = imfilter( curvature(open), Kernel.gaussian((σκ,)) )
    cutline1 === nothing && return cutline2
    cutline2 === nothing && return cutline1
    # indices into open curve at one edge of cut
    open_i1, open_i2 = cutline1[3], cutline2[3]
    # select the cut with lower curvature in the open curve
    #@show (κ[open_i1], κ[open_i2])
    κ[open_i1] < κ[open_i2] ? cutline1 : cutline2
end

function adjust_cutline( cutline, outer, inner, open, d=5, σκ=2 )
    κout = imfilter( curvature(outer), Kernel.gaussian((σκ,)), "circular" )
    κin = imfilter( curvature(inner), Kernel.gaussian((σκ,)), "circular" )
    nout, nin = length(outer), length(inner)

    i_out, i_in, _ = cutline
    out_range, in_range = i_out-d:i_out+d, i_in-d:i_in+d
    i_out = argmin(Dict(i=>κout[mod1(i,nout)] for i in out_range))
    i_in = argmax(Dict(i=>κin[mod1(i,nin)] for i in in_range))

    i_out, i_in
end

function joincurves( outer, inner, i_out, i_in )
    # Assume both curves are positively oriented
    Closed2DCurve([outer.vertices[1:i_out]; inner.vertices[i_in:-1:1];
                   inner.vertices[end:-1:i_in+1]; outer.vertices[i_out+1:end]])
end


## Plotting curves

using RecipesBase
using GeometryBasics: Point

@recipe f(c::Closed2DCurve) = vertices_closed(c)
