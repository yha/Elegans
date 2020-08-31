using OffsetArrays
using GeometryBasics

# TODO: missingpoints can be changed from NaN to missing when ImageFiltering
#       supports missing values.
const missingpoint = Point(NaN, NaN)
ismissingpoint(p) = all(isnan,p)


# trajectories of worm ends end center in "global" coordinates
# (coordinates used in `traj`).
function end_trajectories( traj, contours, irange, midframe )
    # split contours in frame coordinates
    splits_fr, ratios = contour_cache_aligned_split(contours, irange)
    ends_found = .!ismissing.(ratios)
    # n_contours = OffsetArray(
    #     [c isa AbstractVector ? length(c) : 0 for c in contours.(irange)],
    #     irange)

    splits_fr_ok = splits_fr[ends_found]

    miss2nan(x) = replace(x, missing=>NaN)
    x, y = miss2nan( traj.x[irange]), miss2nan(traj.y[irange] )
    center = OffsetArray( Point.(x, y), irange )
    shift = center .- [midframe]
    shift_ok = shift[ends_found]

    # split contours in "global" coordinates
    #@show shift_ok[begin] splits_fr_ok[begin]
    #@show axes(shift_ok,1) axes(splits_fr_ok,1)
    splits_ok = [(c1 .+ s, c2 .+ s) for ((c1,c2),s) in zip(splits_fr_ok, shift_ok)]

    splits = spread( splits_ok, ends_found )

    # worm ends
    ends_ok = [s[1][[1, end]] for s in splits_ok]
    e_ok1 = [first(s[1]) for s in splits_ok]
    e_ok2 = [last(s[1]) for s in splits_ok]

    e1 = spread( e_ok1, ends_found, missingpoint )
    e2 = spread( e_ok2, ends_found, missingpoint )
    #d = e2 - e1

    # @progress "splines" splines_ok =
    #     [try_return(() -> line2spline.(spl)) for spl in splits_ok]
    # splines = Elegans.spread(splines_ok, ends_found)


    # e1 = e_fr1 .+ shift
    # e2 = e_fr2 .+ shift

    e1, e2, center, splits, ratios#, splines
end

##

using IterTools
using StatsBase
using ImageFiltering
using ImageFiltering.KernelFactors: gaussian
using PolygonOps, LinearAlgebra
using OffsetArrays: no_offset_view

# given bool array `v`, return bool array which is true inside sequences of
# `true` of length at least `maxskip+1` in `v`.
# This requires ImageFiltering version ≥ 0.6.14 to work correctly on offset arrays
insequence(v, maxskip) = mapwindow(any, mapwindow(all, v, 0:maxskip, border=Fill(false)), -maxskip:0)

area(c::Closed2DCurve) = PolygonOps.area(vertices_closed(c))
arclen(c::Closed2DCurve) = sum(norm.(diff(vertices_closed(c))))
roundness(c) = 4π*area(c) / arclen(c)^2

function end_assignment_segments( ratios, contours, max_skip=1, max_ratio=0.2,
                                  max_round_z=3, round_winlen=10_001 )
    irange = axes(ratios,1)

    rnd = OffsetArray( map( c -> c isa AbstractVector && !isempty(c) ? roundness(c[1]) : missing,
                       contours.(irange) ), irange )
    #mw(f) = mapwindow( f∘skipmissing, rnd, round_winlen )
    @info "computing roundness statistics"
    @time meanstd_rnd = mapwindow(rnd, round_winlen) do w
        sw = skipmissing(w)
        μ = mean(sw)
        σ = std(sw, mean=μ)
        (μ = μ, σ = σ)
    end
    rnd_μ = (x.μ for x in meanstd_rnd)
    rnd_σ = (x.σ for x in meanstd_rnd)
    rnd_z = (rnd .- rnd_μ) ./ rnd_σ

    too_far = coalesce.(ratios, 0) .> max_ratio
    too_round = coalesce.(rnd_z, 0) .> max_round_z
    too_many_missing_frames = insequence(ismissing.(ratios), max_skip)
    n_contours = OffsetArray(
        [c isa AbstractVector ? length(c) : 0 for c in contours.(irange)],
        irange)

    frames_to_skip = too_round .| too_far .| too_many_missing_frames .| (n_contours .> 1)

    vals, lens = rle(no_offset_view(frames_to_skip))
    # add anything before irange to the first skip range
    if vals[1]
        lens[1] += first(irange)-1
    else
        pushfirst!(vals,true)
        pushfirst!(lens,first(irange)-1)
    end
    clens = cumsum(lens)
    ranges = [range(from + 1, stop = to)
              for (from, to) in partition(clens[1+!vals[1]:end], 2)]
    skip_ranges = [range(from + 1, stop = to)
              for (from, to) in partition(clens[1+vals[1]:end], 2)]

    ranges, skip_ranges, frames_to_skip, rnd, rnd_μ, rnd_σ
end
