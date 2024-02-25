using OffsetArrays
using GeometryBasics

# TODO: missingpoints can be changed from NaN to missing when ImageFiltering
#       supports missing values.
const missingpoint = Point(NaN, NaN)
ismissingpoint(p) = all(isnan,p)


# trajectories of worm ends and center in "global" coordinates
# (coordinates used in `traj`).
# TODO exceptions recorded by `contour_cache_aligned_split` are not preserved
# here (replaced by `missing`)
function end_trajectories( traj, contours, irange )
    # split contours in frame coordinates
    splits_fr, ratios = contour_cache_aligned_split(contours, irange)
    ends_found = .!ismissing.(ratios)

    splits_fr_ok = splits_fr[ends_found]

    miss2nan(x) = replace(x, missing=>NaN)
    x, y = miss2nan( traj.x[irange]), miss2nan(traj.y[irange] )
    center = OffsetArray( Point.(x, y), irange )
    center_ok = center[ends_found]

    # split contours in "global" coordinates
    splits_ok = [(c1 .+ p, c2 .+ p) for ((c1,c2),p) in zip(splits_fr_ok, center_ok)]

    splits = spread( splits_ok, ends_found )

    # worm ends
    ends_ok = [s[1][[1, end]] for s in splits_ok]
    e_ok1 = [first(s[1]) for s in splits_ok]
    e_ok2 = [last(s[1]) for s in splits_ok]

    e1 = spread( e_ok1, ends_found, missingpoint )
    e2 = spread( e_ok2, ends_found, missingpoint )

    e1, e2, center, splits, ratios
end
