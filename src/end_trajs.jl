using OffsetArrays
using GeometryBasics

# TODO: missingpoints can be changed from NaN to missing when ImageFiltering
#       supports missing values.
const missingpoint = Point(NaN, NaN)
ismissingpoint(p) = all(isnan,p)


# trajectories of worm ends end center in "global" coordinates
# (coordinates used in `traj`).
function end_trajectories( traj, contours, irange )
    # split contours in frame coordinates
    splits_fr, ratios = contour_cache_aligned_split(contours, irange)
    ends_found = .!ismissing.(ratios)
    # n_contours = OffsetArray(
    #     [c isa AbstractVector ? length(c) : 0 for c in contours.(irange)],
    #     irange)

    splits_fr_ok = splits_fr[ends_found]

    miss2nan(x) = replace(x, missing=>NaN)
    x, y = miss2nan( traj.x[irange]), miss2nan(traj.y[irange] )
    center = OffsetArray( Point.(x, y), irange )
    center_ok = center[ends_found]

    # split contours in "global" coordinates
    #@show shift_ok[begin] splits_fr_ok[begin]
    #@show axes(shift_ok,1) axes(splits_fr_ok,1)
    splits_ok = [(c1 .+ p, c2 .+ p) for ((c1,c2),p) in zip(splits_fr_ok, center_ok)]

    splits = spread( splits_ok, ends_found )

    # worm ends
    ends_ok = [s[1][[1, end]] for s in splits_ok]
    e_ok1 = [first(s[1]) for s in splits_ok]
    e_ok2 = [last(s[1]) for s in splits_ok]

    e1 = spread( e_ok1, ends_found, missingpoint )
    e2 = spread( e_ok2, ends_found, missingpoint )
    #d = e2 - e1

    # @progress "splines" splines_ok =
    #     [try_return(() -> line2spline.(spl)) for spl in splits_ok]
    # splines = Elegans.spread(splines_ok, ends_found)


    # e1 = e_fr1 .+ shift
    # e2 = e_fr2 .+ shift

    e1, e2, center, splits, ratios#, splines
end
