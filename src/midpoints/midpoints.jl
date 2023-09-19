using OffsetArrays
using Missings
using ImageFiltering
using LinearAlgebra
using Unzip

_mid(sp,s) = Point2(Elegans.splines2midpoint(sp[1],sp[2],s))

"""
Compute midpoints in a range of frames `irange`.
This computation cannot be performed on each frame separately, as the head and tail are identified 
based on their movement. The range should typically be an entire stage.
Returns a named tuple of offset arrays `(; midpts, iters, conf)` where `midpts[i,frame]` is the 
midpoint at frame `frame` and `i`th node (`s=s[i]`), `conf[frame]` is the confidence in the assignment of 
head and tail, and `iters[frame]` is the number of iteration until spline resampling converged.
"""
function range_midpoints( traj, contours, irange, s=0:0.025:1,
                            headtail_method=SpeedHTCM(5,0), end_assigment_params=EndAssigmentParams();
                            resample_max_iter = 100
                        )
    @info "Locating head and tail..."
    head, tail, splits_ht, conf = headtail_trajectories( traj, contours, irange, headtail_method, end_assigment_params )
    @info "Creating contour splines..."
    @time @progress "Fitting splines..." splines = [passmissing(line2spline).(split) for split in splits_ht]

    @info "Computing midlines..."
    @time @progress "midlines" midpts = [passmissing(_mid)(splines[i],t) for i in irange, t in s]

    # TODO compute directly with Point(NaN,NaN) rather than `missing` in previous stages
    midpts = replace(midpts, missing => Elegans.missingpoint)

    # TODO resampling each contour spline first would probably be a bit better 
    # than just resampling the midline
    @info "Resampling mid-points..."
    max_iters = resample_max_iter
    @time midpt_vecs, iters = unzip(Elegans.resample_line(m, s; max_iters, warn_not_converged = false) 
                                     for m in eachrow(midpts))
    not_converged = findall(>(max_iters), iters)
    if !isempty(not_converged)
        @warn "Resampling failed to converge on frames: $(irange[not_converged])" 
    end

    midpts = permutedims(reduce(hcat, midpt_vecs))

   (; midpts = OffsetArray(midpts, irange, eachindex(s)), 
      iters = OffsetArray(iters, irange),
      conf
   )
end

# used instead of skipmissing as a workaround for Statistics.jl issue #50
nomiss(v) = disallowmissing(filter(!ismissing,v))
_cov(x::AbstractVector{P}) where {N,T,P<:Point{N,T}} = isempty(x) ? fill(T(NaN),SMatrix{N,N}) : cov(x)

function make_windows(irange, winlen=60, dwin=20)
    fi = firstindex(irange)
    [irange[i0:i0+winlen-1] for i0 in fi:dwin:length(irange)+fi-winlen]
end

function midpoint_covs( midpts, t, windows )
    covs = [_cov(nomiss(midpts[i,j])) for i in windows, j in eachindex(t)]
    covs, t
end

function midpoint_covs_for_stage( ex, cam, mids, stage_i; winlen=60, dwin=20, stagedict=loadstages() )
    irange = stage_frames( ex, cam, stage_i; stagedict )
    @info "$cam, stage $stage_i: frames $irange"
    midpts = mids(irange)

    t = mids.t

    covs, windows, t = midpts_covs( midpts, t, irange; winlen, dwin )

    (; covs, windows, t, irange)
end

curvelen(pts) = sum(norm(pts[i+1]-pts[i]) for i in firstindex(pts):lastindex(pts)-1)

win_mean_curvelens( pts, windows ) = [mean(curvelen.(eachrow(view(pts,w,:)))) for w in windows]

function normed_midpoint_covs( midpts, t, windows )
    covs, _ = midpoint_covs( midpts, t, windows )
    mean_lens = win_mean_curvelens( midpts, windows )
    normed_covs = covs ./ coalesce.(mean_lens,NaN)
end

const newaxis = [CartesianIndex()]
const newaxis_c = OffsetVector([CartesianIndex()], 0:0)
const D = centered([-1 / 2, 0, 1 / 2]) # derivative
const rD = .- D # reverse derivative


_m2mp(a::AbstractArray{<:Union{Point,Missing}}) = replace(a, missing => Elegans.missingpoint)
_m2mp(a::AbstractArray{<:Point}) = a

# velocity vector at each midpoint. 
# `midpts[i,j]::Point2` is the midpoint in frame i at position s[j] along the midline.
function velocities(midpts)
    midpts_nan = _m2mp(midpts)
    imfilter(midpts_nan, D[:, newaxis_c])
end

# length of midline at each frame.
# `midpts[i,j]::Point2` is the midpoint in frame i at position s[j] along the midline.
function midline_lengths(midpts)
    midpts_nan = _m2mp(midpts)
    curvelen.(eachrow(midpts_nan))
end

# filter to apply to a vector `v` of a positions at consecutive frames 
# to obtain its derivative (velocity) after filtering with `prefilter`:
# (v ⋆ prefilter) ⋆ D = v ⋆ (prefilter ⋆ reverse(D))
# where "⋆" is cross-correlation (imfilter)
# FIXME: add the required zero-padding for this to workaround
velocities_filter(prefilter) = imfilter(prefilter, rD)
velocities_filter(::Nothing) = D

# velocities filter operating per row
velocitives_row_filter(prefilter) = velocities_filter(prefilter)[:, newaxis_c]

# velocity vector at each midpoint. 
# `midpts[i,j]::Point2` is the midpoint in frame i at position s[j] along the midline.
function velocities2(midpts, prefilter = nothing)
    midpts_nan = _m2mp(midpts)
    filt = velocitives_row_filter(prefilter)
    imfilter(midpts_nan, filt)
end

function normal_speeds(midpts; v = nothing)
    midpts_nan = _m2mp(midpts)
    # tangents pointing towards tail
    tangents = LinearAlgebra.normalize.(imfilter(midpts_nan, D[newaxis_c, :]))
    # normals pointing to the worm's right
    nrm = [[-p[2],p[1]] for p in tangents]

    if v === nothing
        v = velocities(midpts_nan)
    end
    #v_t = dot.(tangents, v)
    v_n = dot.(nrm, v)
end

normed_normal_speeds(midpts; v = nothing) = normal_speeds(midpts; v) ./ coalesce.(curvelen.(eachrow(midpts)), NaN)

# Turtle-graphics representations

# workaround for julia issue #32888
_rem2pi(x, mode) = isfinite(x) ? rem2pi(x, mode) : NaN

function curve2turtle(x, y)
    dx, dy = diff(x), diff(y)
    d = hypot.(dx, dy)
    # Zygote doesn't like comprehensions (issue #804)
    angles = map((x, y) -> atan(y, x), dx, dy)
    turns = _rem2pi.(diff(angles), RoundNearest)
    (; d, turns)
end
function curve2turtle(points::AbstractVector{<:Point2})
    xy = reinterpret(reshape, eltype(eltype(points)), points)
    @views curve2turtle(xy[1, :], xy[2, :])
end

# turtle-graphics at each from.
# `midpts[i,j]::Point2` is the midpoint in frame i at position s[j] along the midline.
mids2turtle(midpts) = mapslices(c->Elegans.curve2turtle(c), midpts; dims=2)
mids2turns(midpts) = mapslices(c->Elegans.curve2turtle(c).turns, midpts; dims=2)
