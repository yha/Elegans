using Elegans
using Images, ImageFiltering
using StaticArrays, OffsetArrays
using LinearAlgebra
using GeometryBasics
using Statistics
using Missings

# TODO some functions here should be moved into the package
function load_cam( root, ex, cam, contours_path, midpoints_path;
                    contour_method=Thresholding(1.0,0.34),
                    headtail_method=Elegans.SpeedHTCM(5,0), end_assignment_params=Elegans.EndAssigmentParams() )
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

# for use when saved midpoints are already available
function cam_ncovs( ex, cam, midpoints_path;
                    contour_method=Thresholding(1.0,0.34), headtail_method=Elegans.SpeedHTCM(5,0),
                    winlen=60, dwin=20, t=0:0.025:1 )
#    @info "Loading midpoints for $ex: $cam"
    midsdict = Elegans.load_midpoints( Elegans.midpoints_filename( ex, cam, t;
                                midpoints_path, contour_method, headtail_method,
                                end_assignment_params=Elegans.EndAssigmentParams() ) )
#    @info "Computing normalized midline covariances ($ex: $cam)"
    Dict(irange => normed_midpoint_covs( midpts, t, Elegans.make_windows(irange, winlen, dwin) ) 
                for (irange,midpts) in midsdict)
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
