using OffsetArrays
using Missings
using ImageFiltering
using LinearAlgebra

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

# used instead of skipmissing as a workaround for Statistics.jl issue #50
nomiss(v) = disallowmissing(filter(!ismissing,v))
_cov(x::AbstractVector{P}) where {N,T,P<:Point{N,T}} = isempty(x) ? fill(T(NaN),SMatrix{2,2}) : cov(x)

function midpoint_covs( midpts, t, irange = axes(midpts,1); winlen=60, dwin=20 )
    fi = firstindex(irange)
    windows = [irange[i0:i0+winlen-1] for i0 in fi:dwin:length(irange)+fi-winlen]
    covs = [_cov(nomiss(midpts[i,j])) for i in windows, j in eachindex(t)]
    covs, windows, t
end

function midpoint_covs_for_stage( ex, cam, mids, stage_i; winlen=60, dwin=20, stagedict=loadstages() )

    irange = stage_frames( ex, cam, stage_i; stagedict )
    @info "$cam, stage $stage_i: frames $irange"
    midpts = mids(irange)

    t = mids.t

    covs, windows, t = midpts_covs( midpts, t, irange; winlen, dwin )

    (;covs, windows, t, irange)
end


function normalize_covs( covs, midpts, windows )
    curvelen(pts) = sum(norm.(diff(pts)))
    lens = curvelen.(eachrow(midpts))
    mean_lens = [mean(lens[w]) for w in windows]

    covs ./ coalesce.(mean_lens,NaN)
end

# function normed_stg_midpts_covs( ex, cam, mids, stage_i; winlen=60, dwin=20 )
#     covs, windows, t, irange = midpts_covs_for_stage( ex, cam, mids, stage_i )
#
#     midpts = mids(irange)
#     normed_covs = normalize_covs( covs, midpts, windows )
#     normed_covs, windows, t, irange
# end
#
function normed_midpoint_covs(midpts, t, irange = axes(midpts,1); winlen=60, dwin=20)
    covs, windows, _ = midpoint_covs( midpts, t, irange; winlen, dwin )
    normed_covs = normalize_covs( covs, midpts, windows )
    normed_covs, windows
end


const newaxis = [CartesianIndex()]
const newaxis_c = OffsetVector([CartesianIndex()],0:0)

function normal_speeds(midpts)
    d = centered([-1/2,0,1/2]) # derivative

    midpts_nan = replace(midpts, missing=>Elegans.missingpoint)
    # tangents pointing towards tail
    tangents = LinearAlgebra.normalize.(imfilter(midpts_nan, d[newaxis_c,:]))
    # normals pointing to the worm's right
    nrm = [[-p[2],p[1]] for p in tangents]

    v = imfilter(midpts_nan, d[:,newaxis_c])
    #v_t = dot.(tangents, v)
    v_n = dot.(nrm, v)
end
