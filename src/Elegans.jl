module Elegans

export import_coords, import_and_calc, loadstages, stage_names, stage_frames,
       VideoCache, get_frame, nframes, video_index,
       mark_stages_gui, mark_stages_window, frames_per_s,
       raw_worm_contours, worm_contour, curvature,
       Thresholding, SeededSegmentation,
       incurve, cutline, joincurves,
       Closed2DCurve, n_highest_peaks_circular,
       contour_cache, init_contours, save_contours,
       midline_by_thinning, line2spline, resample_spline, spline_anglevec,
       anglevec_by_thinning, midline_cache, spline_cache, points,
       ends_alignment_mask, align_ends!, find_end_indices,
       contour2splines, aligned_splines,
       isroaming, roam_for_stage,
       end_trajectories, end_assignment_segments, headtail_trajectories,
       range_midpoints, midpoint_cache, init_midpoints,
       save_midpoints, midpoints_filename,
       forward_speed, speed_heatmap, speed_heatmap_data, speed_heatmap_plot,
       log_speed_ratios,
       speed_stats, direction_stats, spread_stats,
       plot_frame_with_time,
       # TODO move to a new packge?
       try_return, passex, missex

include("utils.jl")
include("mark_stages.jl")
include("contour.jl")
include("segment_contour.jl")
include("midline.jl")
include("shape.jl")
include("peaks.jl")
include("videocache.jl")
include("caching.jl")
include("end_trajs.jl")
include("roam_fraction.jl")
include("video_output.jl")
include("headtail/forward.jl")
include("headtail/lsr.jl")
include("headtail/traj_stats.jl")
include("headtail/headtail_traj.jl")


end # module
