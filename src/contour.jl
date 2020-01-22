using Contour
import Images, ImageFiltering
using Images: otsu_threshold, imfilter
using CircularArrays
#using OffsetArrays
using GeometryTypes
using LinearAlgebra
import Luxor
import IterTools
import ShiftedArrays
const cshift = ShiftedArrays.circshift

# A closed 2d curve with circular indexing
struct Closed2DCurve{T} <: AbstractVector{Point{2,T}}
    vertices::CircularVector{Point{2,T}}
    function Closed2DCurve(v::AbstractVector{Point{2,T}}; orient_pos=false) where T
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

function Closed2DCurve( c::Contour.Curve2; orient_pos=false )
    _is_closed(c) || throw(ArgumentError("Curve not closed"))
    # Contour.Curve2 repeats the first vertex when closed.
    # Remove the repeated vertex.
    Closed2DCurve( Point.(c.vertices[1:end-1]); orient_pos=orient_pos )
end
vertices(c::Closed2DCurve) = c.vertices
vertices_closed(c::Closed2DCurve) = push!(Vector(c.vertices),first(c.vertices))

Base.size(c::Closed2DCurve) = (length(c.vertices),)
Base.getindex(c::Closed2DCurve, i::Int) = c.vertices[i]

Base.iterate(c::Closed2DCurve, state=firstindex(c.vertices)) = ( state > lastindex(c.vertices)
                ? nothing : (c.vertices[state], state+1) )
#Base.eltype(::Type{Closed2DCurve{T}}) where T = Point{2,T}
#Base.IteratorSize(::Type{Closed2DCurve{T}}) where T = Base.HasLength()

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

#angle(p::Point2) = atan(p[2],p[1])

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

# ##
#
# function _circshift(c::Curve2, shift)
#     # remove the last vertex which repeats the first
#     v = c.vertices[1:end-1]
#     vs = circshift(v,shift)
#     # close the curve by adding the first vertex as the last
#     push!(vs,vs[1])
# end
# function Base.circshift!(c::Curve2, shift)
#     c.vertices .= _circshift(c,shift)
#     c
# end
# Base.circshift(c::Curve2, shift) = Curve2(_circshift(c,shift))
# Base.reverse!(c::Curve2) = reverse!(c.vertices)
#
# function _circfilter(c::Curve2, kernel)
#     v = imfilter( c.vertices[1:end-1], kernel, "circular" )
#     push!(v,v[1])
# end
# circfilter(c::Curve2, kernel) = Curve2(_circfilter(c,kernel))
# function circfilter!(c::Curve2, kernel)
#     c.vertices .= _circfilter(c,kernel)
#     c
# end

_is_closed(c::Curve2) = c.vertices[end] == c.vertices[1]

function raw_worm_contours(img::AbstractArray{<:AbstractFloat}, level::Number)::Vector{Closed2DCurve{Float64}}
    # Contour.contour(x,y,z) expects z to be indexed as
    # z[x[i],y[i]], contrary to the Images.jl convention img[y,x].
    # To get coordinates ordered as (x,y), we transpose img
    y, x = float.(axes(img))
    curves = lines(Contour.contour(x, y, img', level))
    curves = filter(_is_closed, curves)
    # # Remove last repeated point (same as first since we picked only closed curves),
    # # as well as consecutive duplicates, from each curve
    # curves = [Point.(StatsBase.rle(curve.vertices)[1][1:end-1]) for curve in curves]
    # filter!( !isempty, curves )
    #filter!( c->!isempty(c.vertices), curves )
    closed_curves = Closed2DCurve.(curves; orient_pos=true)
    sort( closed_curves; by=length, rev=true )
end
raw_worm_contours(img::AbstractArray{<:Number}, level::Number) = raw_worm_contours(float.(img), level)
raw_worm_contours(img::AbstractArray{<:Gray}, level) = raw_worm_contours(real.(img), real(level))

# function raw_worm_contour(img, level)
#     curves = raw_worm_contours(img,level)
#     #length(curves) == 1 || @warn("Contour has multiple components")
#     curves[1]
# end

# curvature(c::Curve2) = curvature(coordinates(c)...)
# function _curve_stats(x,y)
#     (x[1],y[1]) == (x[end],y[end]) || error("Curve not closed")
#     dx,dy = diff(x), diff(y)
#     dl = hypot.(dx,dy)
#     dl_mid = (dl .+ circshift(dl,1)) ./ 2
#     θ = atan.(dy,dx)
#     θ = [θ[end]; θ]
#     dθ = fix_angle.(diff(θ))
#     dl_mid, dθ
# end
# function curvature(x,y)
#     dl_mid, dθ = _curve_stats(x,y)
#     if any(dl_mid .< 1e-3)
#         @warn("Some curve points are very close ($(minimum(dl_mid)) pixels " *
#               "between mid-segments along curve)", maxlog=4)
#     end
#     κ = dθ ./ dl_mid # curvature
# end

# function worm_contour(img, level=otsu_threshold(img))
#     c = raw_worm_contour(img, level)
#     # place the highest (positive signed) curvature element first
#     κ = curvature(c)
#     i_max = argmax(κ)
#     circshift!(c, 1-i_max)
# end

# Curve predicates, using Luxor.jl

luxorpoint(p::Luxor.Point) = p
luxorpoint(p) = Luxor.Point(p[1],p[2])

function incurve( p, c::Closed2DCurve )
    Luxor.isinside( luxorpoint(p), luxorpoint.(c); allowonedge=true )
end

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
using GeometryTypes: Point

# @recipe f(c::Curve2) = Point.(c.vertices)
@recipe f(c::Closed2DCurve) = vertices_closed(c)
#@recipe f(::Type{<:Closed2DCurve}, ::Closed2DCurve) = (vertices_closed, identity)
