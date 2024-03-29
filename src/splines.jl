using Dierckx
using IterTools

splines2midline(spl1, spl2) = x -> splines2midpoint.(Ref(spl1), Ref(spl2), x)
splines2midpoint(spl1, spl2, x) = mean((spl1(x),spl2(x)))

merge_nearby_points(line, th=0.01) = IterTools.imap( mean,
                                group_pairwise( (x,y) -> norm(x-y)<th, line ) )

function line2spline(line, th=0.01)
    # Merge almost-repeated points (which are sometimes generated by the contour
    # computation), to avoid numerical issues in spline computation.
    line = collect(merge_nearby_points(line,th))
    d = norm.(Tuple.(diff(line)))
    d ./= sum(d)
    # length parameter
    s = [0.0; cumsum(d)]
    # # vector of CartesianIndex => array with [x,y] pairs in columns
    # A = mapreduce(i->[i[2],i[1]], hcat, line)
    #A = reduce(hcat, line)
    A = reinterpret(reshape, eltype(eltype(line)), line)

    # output spline is (very) approximately length-parameterized
    spl = ParametricSpline(s, A; k = min(3,length(s)-1))
end

resample_spline(spl, s) = line2spline(eachcol(spl(s)))


points(spline, x) = Point2.(eachcol(spline(x)))

resample_line_once(line, s) =
    any(p -> any(!isfinite, p), line) ?
        fill(Elegans.missingpoint, size(line)) :
        points( line2spline(line), s )

# Iteratively fit spline and resample until distance of adjacent points is nearly constant
# Returns the `(line, iters)` where `line` is the resampled line and `iters` is the 
# number of iterations until convergence, or max_iters+1 if not converged
function resample_line(line, s; threshold = 0.99, max_iters = 100, warn_not_converged = true)
    evenness(x) = let (d_min, d_max) = extrema(norm, diff(x))
        d_min / d_max
    end
    iters = 0
    while evenness(line) < threshold
        iters += 1
        if iters > max_iters
            warn_not_converged && @warn "Spline resampling failed to converge after $max_iters iterations." evenness=evenness(line)
            break
        end
        line = resample_line_once(line, s)
    end
    line, iters
end

