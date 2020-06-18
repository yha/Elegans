# function mean_var_nonan(e)
#     ef = filter(p->!any(isnan,p),e)
#     μ, Σ = mean(ef), cov(ef)
# end
#
# function dir_circvar(e)
#     dirs = 2 .* angle.(diff(v))
#     circvar.(filter.(isfinite,dirs))
# end

## speed

using ImageFiltering: centered

cdiff(v) = imfilter(v, centered([-1/2,0,1/2]))#, Inner())

function speed_stats(e1, e2, center, t, s)
    k = gaussian(t)
    n1, n2, n_c = (norm.(cdiff(imfilter(v,k,NA()))) for v in (e1, e2, center))
    lsr = imfilter(log.(n1 ./ n2), gaussian(s), NA())
    lsr, n1, n2, n_c
end


#using StaticArrays

_cov(v::AbstractArray{<:StaticVector{N}}) where N = isempty(v) ? (@SMatrix fill(NaN,N,N)) : cov(v)
cov_finite(v) = _cov(filter(p -> any(isfinite, p), v))
function mean_cov_finite(v)
    vf = filter(p -> any(isfinite, p), v)
    (mean = mean(vf), cov = _cov(v))
end

bounded_log_ratio( x1, x2, lbound = 1e-12) = log( max(x1,lbound) / max(x2,lbound) )

function spread_stats(e1, e2, center, winlen, s, mindet=1e-12)
    d1, d2, d_c = (mapwindow( det ∘ cov_finite, v, winlen ) for v in (e1, e2, center))
    lvr = imfilter( bounded_log_ratio.( d1, d2, mindet ), gaussian(s), NA() )
    lvr, d1, d2, d_c
end



## directional variance

Base.angle(p::Point2) = atan(p[2],p[1])

#circmean(v) = angle(mean(exp.(v.*im)))
#circvar(v) = 1 - norm(mean(exp.(v.*im)))
function direction_stats(e1, e2, center, t, winlen, s)
    k = gaussian(t)
    a(v) = angle.(cdiff(imfilter(v,k,NA())))
    a1 = a(e1)
    a2 = a(e2)
    a_c = a(center)
    dirvar(w) = circvar(2 .* w)
    c1 = mapwindow( dirvar, a1, winlen )
    c2 = mapwindow( dirvar, a2, winlen )
    c_c = mapwindow( dirvar, a_c, winlen )
    ldvr = imfilter( log.( c1 ./ c2 ), gaussian(s), NA() )
    ldvr, c1, c2, c_c, a1, a2, a_c
end

# polar(p::Point2) = (norm = norm(p), angle=angle(p))
# from_polar(polar) = from_polar( polar.norm, polar.angle )
# from_polar(norm, angle) = Point2( norm*cos(angle), norm*sin(angle) )
#
# function direction_stats1(e1, e2, center, t, winlen)
#     k = gaussian(t)
#     d(v) = polar.(diff(imfilter(v,k,NA())))
#     d1 = d(e1)
#     d2 = d(e2)
#     d_c = d(center)
#     dirvar_weighted(w) = 1 ./ norm(mean( from_polar(n,2a) for (n,a) in w ))
#     #dirvar_weighted(w) = 1 ./ norm(mean( from_polar(1,2a) for (n,a) in w ))
#     c1 = mapwindow( dirvar_weighted, d1, winlen )
#     c2 = mapwindow( dirvar_weighted, d2, winlen )
#     c_c = mapwindow( dirvar_weighted, d_c, winlen )
#     ldvr = bounded_log_ratio.( c1, c2 )
#     ldvr, c1, c2, c_c
# end
