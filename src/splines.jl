using Dierckx
using IterTools

merge_nearby_points(line, th=0.01) = IterTools.imap( mean,
                                group_pairwise( (x,y) -> norm(x-y)<th, line ) )

function line2spline(line)
    line = merge_nearby_points(line)
    d = norm.(Tuple.(diff(line)))
    d ./= sum(d)
    # length parameter
    s = [0.0; cumsum(d)]
    # # vector of CartesianIndex => array with [x,y] pairs in columns
    # A = mapreduce(i->[i[2],i[1]], hcat, line)
    A = reduce(hcat, line)

    # output spline is approximately length-parameterized (exact at knot points)
    spl = ParametricSpline(s,A)
end

resample_line(line, s) = Point2.(eachcol(line2spline(line)(s)))

resample_spline(spl, s) = ParametricSpline(s,spl(s))
