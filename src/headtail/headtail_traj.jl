using IterTools
using StatsBase
using ImageFiltering
using ImageFiltering.KernelFactors: gaussian
using GeometryBasics
using PolygonOps
using LinearAlgebra
using OffsetArrays
using OffsetArrays: no_offset_view
using RLEVectors
using Missings


Base.@kwdef struct EndAssigmentParams
    max_skip::Int = 1
    max_ratio::Float64 = 0.2
    max_round_z::Float64 = 3
    round_winlen::Int = 10_001
end

# given bool array `v`, return bool array which is true inside sequences of
# `true` of length at least `maxskip+1` in `v`.
# This requires ImageFiltering version ≥ 0.6.14 to work correctly on offset arrays
insequence(v, maxskip) = mapwindow(any, mapwindow(all, v, 0:maxskip, border=Fill(false)), -maxskip:0)

area(c::Closed2DCurve) = PolygonOps.area(vertices_closed(c))
arclen(c::Closed2DCurve) = sum(norm.(diff(vertices_closed(c))))
roundness(c) = 4π*area(c) / arclen(c)^2

function end_assignment_segments( ratios, contours, params=EndAssigmentParams() )
    irange = axes(ratios,1)

    rnd = OffsetArray( map( c -> c isa AbstractVector && !isempty(c) ? roundness(c[1]) : missing,
                       contours.(irange) ), irange )
    #mw(f) = mapwindow( f∘skipmissing, rnd, round_winlen )
    @info "computing roundness statistics"
    @time meanstd_rnd = mapwindow(rnd, params.round_winlen) do w
        sw = skipmissing(w)
        μ = mean(sw)
        σ = std(sw, mean=μ)
        (μ = μ, σ = σ)
    end
    rnd_μ = (x.μ for x in meanstd_rnd)
    rnd_σ = (x.σ for x in meanstd_rnd)
    rnd_z = (rnd .- rnd_μ) ./ rnd_σ

    too_far = coalesce.(ratios, 0) .> params.max_ratio
    too_round = coalesce.(rnd_z, 0) .> params.max_round_z
    too_many_missing_frames = insequence(ismissing.(ratios), params.max_skip)
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

function headtail_splits( e1, e2, splits, y )
    axes(e1) == axes(e2) == axes(splits) == axes(y) || throw(
                ArgumentError("Incompatible axes: $(axes.((e1,e2,splits,y)))"))
    head = similar(e1)
    tail = similar(e1)
    splits_out = similar(splits)
    for i in eachindex(e1)
        if ismissing(splits[i]) || ismissing(y[i])
            head[i] = missingpoint
            tail[i] = missingpoint
            splits_out[i] = missing
        else
            headfirst = y[i] > 0
            head[i] = headfirst ? e1[i] : e2[i]
            tail[i] = headfirst ? e2[i] : e1[i]
            splits_out[i] = headfirst ? splits[i] : reverse.(splits[i])
        end
    end
    head, tail, splits_out, abs.(y)
end

function rlevec(ranges, irange, y)
    firstindex(y) == 1 || throw(ArgumentError("`y` should be 1-indexed"))
    n = length(ranges)
    n == length(y) || throw(DimensionMismatch("`length(y) != length(ranges)`"))

    i0 = first(irange)
    ends_outside_ranges = last(last(ranges)) < last(irange)
    outlen = 2length(ranges) + ends_outside_ranges

    runends = Vector{Int}(undef, outlen)
    runvalues = Vector{Union{eltype(y), Missing}}(undef, outlen)
    for (i,r) in enumerate(ranges)
        runends[2i-1] = first(r) - i0
        runends[2i] = last(r) + 1 - i0
        runvalues[2i-1] = missing  # requires RLEVectors > v0.9.0
        runvalues[2i] = y[i]
    end
    if ends_outside_ranges
        runends[end] = length(irange)
        runvalues[end] = missing
    end

    OffsetVector( RLEVector(runvalues, runends), irange )
end

# trajectories of worm ends end center in "global" coordinates
# (coordinates used in `traj`).
function headtail_trajectories( traj, contours, irange,
                headtail_method=SpeedHTCM(5,0), end_assigment_params=EndAssigmentParams() )
    # Ends (and split contours) aligned by distance ratios
    e1, e2, center, splits_12, ratios = end_trajectories(traj, contours, irange)
    ends_found = .!ismissing.(ratios)

    # Break to segments
    ranges, = Elegans.end_assignment_segments( ratios, contours, end_assigment_params )

    # Head-tail classification
    h, = Elegans.headtail_stats( headtail_method, e1, e2, center )
    mh_per_range = [mean(h[r]) for r in ranges]
    mh = rlevec( ranges, irange, mh_per_range )

    head, tail, splits_ht, conf = Elegans.headtail_splits( e1, e2, splits_12, mh )
end

_mid(sp,s) = Point2(Elegans.splines2midpoint(sp[1],sp[2],s))

function range_midpoints( traj, contours, irange, t=0:0.025:1,
                            headtail_method=SpeedHTCM(5,0), end_assigment_params=EndAssigmentParams() )
    @info "Locating head and tail..."
    head, tail, splits_ht, conf = headtail_trajectories( traj, contours, irange, headtail_method, end_assigment_params )
    @info "Creating contour splines..."
    # # TODO this @progress reqires Juno.jl unmerged PR #605
    # @progress "splines" splines = [passmissing(line2spline).(spl) for spl in splits_ht]
    # firstindex(splines) == 1 && error("Juno.@progress lost array offsets") # Juno should be updated or @progress removed
    # alternative:
    @time splines = [passmissing(line2spline).(spl) for spl in splits_ht]

    @info "Find midpoints..."
    @progress "midpoints" midpts = [passmissing(_mid)(splines[i],t) for i in irange, t in t]
    midpts = OffsetArray(midpts, irange, eachindex(t))
end
