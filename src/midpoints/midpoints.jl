using OffsetArrays
using Missings
using ImageFiltering
using LinearAlgebra
using Unzip

_mid(sp,s) = Point2(Elegans.splines2midpoint(sp[1],sp[2],s))

function range_midpoints( traj, contours, irange, s=0:0.025:1,
                            headtail_method=SpeedHTCM(5,0), end_assigment_params=EndAssigmentParams() )
    @info "Locating head and tail..."
    head, tail, splits_ht, conf = headtail_trajectories( traj, contours, irange, headtail_method, end_assigment_params )
    @info "Creating contour splines..."
    # # TODO this @progress reqires Juno.jl unmerged PR #605
    # @progress "splines" splines = [passmissing(line2spline).(spl) for spl in splits_ht]
    # firstindex(splines) == 1 && error("Juno.@progress lost array offsets") # Juno should be updated or @progress removed
    # alternative:
    @time splines = [passmissing(line2spline).(split) for split in splits_ht]

    @info "Find midpoints..."
    @progress "midpoints" midpts = [passmissing(_mid)(splines[i],t) for i in irange, t in s]

    # TODO compute directly with Point(NaN,NaN) rather than `missing` in previous stages
    midpts = replace(midpts, missing => Elegans.missingpoint)

    @info "Resampling mid-points..."
    max_iters = 100
    @time midpt_vecs, iters = unzip(Elegans.resample_line(m, s; max_iters) for m in eachrow(midpts))
    not_converged = findall(==(max_iters), iters)
    if !isempty(not_converged)
        @warn "Resampling failed to converge on frames: $(irange[not_converged])" 
    end

    midpts = permutedims(reduce(hcat, midpt_vecs))

    OffsetArray(midpts, irange, eachindex(s)), OffsetArray(iters, irange)
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


_m2mp(a::AbstractArray{<:Union{Point,Missing}}) = replace(a, missing => Elegans.missingpoint)
_m2mp(a::AbstractArray{<:Point}) = a

# velocity vector at each midpoint. 
# `midpts[i,j]::Point2` is the midpoint in frame i at position s[j] along the midline.
function velocities(midpts)
    midpts_nan = _m2mp(midpts)
    imfilter(midpts_nan, D[:, newaxis_c])
end


function normal_speeds(midpts)
    midpts_nan = _m2mp(midpts)
    # tangents pointing towards tail
    tangents = LinearAlgebra.normalize.(imfilter(midpts_nan, D[newaxis_c, :]))
    # normals pointing to the worm's right
    nrm = [[-p[2],p[1]] for p in tangents]

    v = velocities(midpts_nan)
    #v_t = dot.(tangents, v)
    v_n = dot.(nrm, v)
end

normed_normal_speeds(midpts) = normal_speeds(midpts) ./ coalesce.(curvelen.(eachrow(midpts)), NaN)