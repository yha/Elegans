using ImageSegmentation
using ImageFiltering
using ImageFiltering.KernelFactors: gaussian
using StatsBase: mode

function find_bg_pixel( frame::AbstractMatrix{<:AbstractGray{<:FixedPoint}} )
    bgcolor = mode(frame)
    findfirst(==(bgcolor),frame)
end

bg_worm_seeded_segment( frame::AbstractMatrix, s=1.0 ) = bg_worm_seeded_segment( Gray.(frame), s )
function bg_worm_seeded_segment( frame::AbstractMatrix{<:AbstractGray}, s=1.0 )
    bg_pixel = find_bg_pixel(frame)
    worm_pixel = argmin(imfilter(frame,gaussian((s,s))))
    seeds = [(worm_pixel,1), (bg_pixel,2)]
    segments = seeded_region_growing(frame, seeds)
end


function bgworm_segment_closed_contours(img, s=1.0)::Vector{Closed2DCurve{Float64}}
    segments = bg_worm_seeded_segment( img, s )
    closed_contour_levels( imfilter(labels_map(segments), gaussian((s,s))), 1.5 )
end
