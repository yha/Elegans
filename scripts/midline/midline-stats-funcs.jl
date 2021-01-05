using Elegans
using Images, ImageFiltering
using StaticArrays, OffsetArrays
using LinearAlgebra
using GeometryBasics
using Statistics
using Missings

##

function load_cam( root, ex, cam, contours_path, midpoints_path;
                    contour_method=Thresholding(1.0,0.34),
                    headtail_method=SpeedHTCM(5,0), end_assignment_params=EndAssigmentParams() )
    relcam = joinpath(ex, cam)
    campath = joinpath(root, relcam)
    @assert isdir(campath)

    contours, contours_file, vcache = init_contours(
                            relcam, root, contour_method, contours_path)

    traj = import_and_calc(relcam, 3, root)

    mids, midfile = Elegans.init_midpoints(ex, cam, traj, contours;
                            contour_method, midpoints_path, headtail_method, end_assignment_params)

    (;traj, contours, mids, vcache, contours_file, midfile)
end

# used instead of skipmissing as a workaround for Statistics.jl issue #50
nomiss(v) = disallowmissing(filter(!ismissing,v))
_cov(x::AbstractVector{P}) where {N,T,P<:Point{N,T}} = isempty(x) ? fill(T(NaN),SMatrix{2,2}) : cov(x)

function framerange( ex, cam, stage_i; stagedict=loadstages() )
    stage_boundaries = stagedict[ex][cam]
    return stage_boundaries[stage_i]+1:stage_boundaries[stage_i+1]
end

function midpts_covs( midpts, t, irange; winlen=60, dwin=20 )
    windows = [irange[i0+1:i0+winlen] for i0 in 0:dwin:length(irange)-winlen]
    covs = [_cov(nomiss(midpts[i,j])) for i in windows, j in eachindex(t)]
    covs, windows, t
end

function midpts_covs_for_stage( ex, cam, mids, stage_i; winlen=60, dwin=20, stagedict=loadstages() )

    irange = framerange( ex, cam, stage_i; stagedict )
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

function normed_stg_midpts_covs( ex, cam, mids, stage_i; winlen=60, dwin=20 )
    covs, windows, t, irange = midpts_covs_for_stage( ex, cam, mids, stage_i )

    midpts = mids(irange)
    normed_covs = normalize_covs( covs, midpts, windows )
    normed_covs, windows, t, irange
end

function normed_midpts_covs(midpts, t, irange; winlen=60, dwin=20)
    covs, windows, _ = midpts_covs( midpts, t, irange; winlen, dwin )
    normed_covs = normalize_covs( covs, midpts, windows )
    normed_covs, windows
end

# for use when saved midpoints are already available
function cam_ncovs( ex, cam, stage_i, midpoints_path;
                    contour_method=Thresholding(1.0,0.34), headtail_method=Elegans.SpeedHTCM(5,0),
                    winlen=60, dwin=20, t=0:0.025:1 )
#    @info "Loading midpoints for $ex: $cam"
    midsdict = Elegans.load_midpoints( Elegans.midpoints_filename( ex, cam, t;
                                midpoints_path, contour_method, headtail_method,
                                end_assignment_params=Elegans.EndAssigmentParams() ) )
    irange = framerange( ex, cam, stage_i )
#    @info "Computing normalized midline covariances ($ex: $cam)"
    Dict((irange, normed_midpts_covs(midpts,t,irange; winlen, dwin)) for (irange,midpts) in midsdict)
end

# for use when saved midpoints may not be available
function cam_ncovs_full( root, ex, cam, stage_i, contours_path, midpoints_path;
                    winlen=60, dwin=20,
                    contour_method=Thresholding(1.0,0.34), headtail_method=Elegans.SpeedHTCM(5,0),
                    stagedict=loadstages() )
    @info "Loading cam $ex: $cam"
    traj, contours, mids, vcache, contours_file, midfile = load_cam(
                                root, ex, cam, contours_path, midpoints_path)
    @info "Computing normalized midline covariances ($ex: $cam)"
    normed_covs, windows, t, irange = normed_midline_covs(ex, cam, mids, stage_i; winlen, dwin)
    @info "Saving midpoints to $midfile"
    save_midpoints(mids, midfile)
    normed_covs, windows, t, irange
end

##

const newaxis = [CartesianIndex()]
const newaxis_c = OffsetVector([CartesianIndex()],0:0)
function normal_speeds(midpts)
    d = centered([-1/2,0,1/2]) # derivative

    midpts_nan = replace(midpts, missing=>Elegans.missingpoint)
    # tangents pointing towards tail
    tangents = normalize.(imfilter(midpts_nan, d[newaxis_c,:]))
    # normals pointing to the worm's right
    nrm = [[-p[2],p[1]] for p in tangents]

    v = imfilter(midpts_nan, d[:,newaxis_c])
    #v_t = dot.(tangents, v)
    v_n = dot.(nrm, v)
end
