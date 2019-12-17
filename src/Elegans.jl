module Elegans

export import_coords, import_and_calc, loadstages, stage_names,
       VideoCache, get_frame, nframes,
       mark_stages_gui, mark_stages_window, frames_per_s,
       raw_worm_counts, worm_contour, curvature, #circfilter, circfilter!,
       incurve, cutline, joincurves,
       Closed2DCurve, n_highest_peaks_circular,
       midline_by_thinning, line2spline, spline_anglevec,
       anglevec_by_thinning, midline_cache, spline_cache,
       ends_alignment_mask, align_ends!

include("mark_stages.jl")
include("shape.jl")
include("contour.jl")
include("midline.jl")
include("peaks.jl")
include("videocache.jl")

end # module
