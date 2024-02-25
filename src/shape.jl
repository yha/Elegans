using Dierckx
using IterTools
using GeometryBasics: Point2
using LinearAlgebra, Statistics
# using ResultTypes
# using ResultTypes: unwrap
using StaticArrays
using OffsetArrays
using ProgressLogging: @progress

include("splines.jl")


function find_end_indices(c::Closed2DCurve, σc = 2.0, σκ = 5.0)
    cf = imfilter(c, Kernel.gaussian((σc,)))
    κ = curvature(cf)

    κm = imfilter(κ, Kernel.gaussian((σκ,)))
    i = collect(keys(n_highest_peaks_circular(κm, 2)))

    if length(i) < 2
        throw(ErrorException("Less than two curvature peaks found"))
    else
        @assert length(i) == 2
        SVector(i...)
    end
end
find_ends( c, σc = 2.0, σκ = 5.0 ) = c[find_end_indices(c,σc,σκ)]

function ends_alignment_mask(ends, return_ratios=false)
    out = Vector{Bool}(undef,length(ends))
    out[1] = false   # Just for consistency. Not necessary for correct output.
    dist(ends1,ends2) = mean(norm.(ends1 .- ends2))
    return_ratios && (ratios = Vector{Float64}(undef,length(ends)-1))
    prev = first(ends)
    for i = firstindex(ends)+1:lastindex(ends)
        curr = ends[i]
        curr_reversed = reverse(curr)
        d, d_reversed = dist(curr, prev), dist(curr_reversed, prev)
        return_ratios && (ratios[i-1] = /(minmax(d,d_reversed)...))
        out[i] = out[i-1] ⊻ (d > d_reversed)
        prev = curr
    end
    if return_ratios
        out, ratios
    else
        out
    end
end

function align_ends!(lines)
    ends = [(l[1], l[end]) for l in lines]
    switch_ends = ends_alignment_mask(ends)
    for i in eachindex(ends)
        switch_ends[i] && (lines[i] = reverse(lines[i]))
    end
    lines
end

function split_curve(c::Closed2DCurve, (i,j) = find_end_indices(c))
    s = i>j
    n = length(c)
    half1, half2 = c[i:j+s*n], c[i+!s*n:-1:j]
end

contour2splines(c::Closed2DCurve, ends=find_end_indices(c)) = line2spline.(split_curve(c,ends))

function swapat!(pairs, mask)
    for p in pairs[mask]
        p[1:2] = [p[2],p[1]]
    end
    pairs
end

function swapat(pairs,mask)
    out = similar(pairs)
    for i in eachindex(pairs)
        out[i] = (mask[i] ? reverse : copy)(pairs[i])
    end
    out
end


"""
Split all contours along aligned ends.
Returns `(splits, ratios)` where
`splits[i]` is either the split for the i-th contour
          or the reported error
`ratios[i]` is the swap ratio for frame i and the next non-failed frame, or
          `missing` on failed frames and on the last non-failed frame
"""
function aligned_split(contours)
    @progress "finding ends" ends_i = [trying(find_end_indices)(c) for c in contours]
    ends = [passex(getindex)(c,i) for (c,i) in zip(contours,ends_i)]
    ok = [!(x isa Exception) for x in ends]
    m, ratios_ok = ends_alignment_mask(ends[ok], true)
    ends_i_aligned = swapat(ends_i[ok],m)
    @progress "splitting curves" splits_ok = [split_curve(c,i)
                            for (c,i) in zip(contours[ok], ends_i_aligned)]

    splits = Vector{Union{eltype(splits_ok),Exception}}(undef,length(contours))
    splits[ok] .= splits_ok
    splits[.!ok] .= ends_i[.!ok]
    splits, spread([ratios_ok;missing],ok)
end

function contour_cache_aligned_split(contours_f, i)
    @progress "fetching contours" cs = [contours_f(i) for i in i]

    contours_found = [x isa Vector && !isempty(x) for x in cs]
    cs1 = [c[1] for c in cs[contours_found]]
    split1, ratios1 = aligned_split(cs1)

    splits = spread(split1,contours_found)
    ratios = spread(ratios1,contours_found)

    OffsetArray(splits,i), OffsetArray(ratios,i)
end


function aligned_splines( frames, contouring_method )
    @progress "contours" cs = [raw_worm_contours(fr, contouring_method) for fr in frames]
    ok = (!isempty).(cs)
    cs1 = [c[1] for c in cs[ok]]
    split1, ratios = aligned_split(cs1)
    @progress "splines" splines = [try_return(()->passex(s->line2spline.(s))(spl)) for spl in split1]
    spread(splines,ok), cs, spread(ratios,ok)
end
