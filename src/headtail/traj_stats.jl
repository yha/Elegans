## speed

using ImageFiltering: centered

cdiff(v) = imfilter(v, centered([-1/2,0,1/2]))#, Inner())

function speed_stats(e1, e2, center, t, s)
    k = gaussian(t)
    n1, n2, n_c = (norm.(cdiff(imfilter(v,k,NA()))) for v in (e1, e2, center))
    lsr = imfilter(log.(n1 ./ n2), gaussian(s), NA())
    lsr, n1, n2, n_c
end


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

##

abstract type HeadTailClassificationMethod end

struct SpeedHTCM{T,S} <: HeadTailClassificationMethod
    t::T
    s::S
end
struct SpreadHTCM{S} <: HeadTailClassificationMethod
    winlen::Int
    s::S
    mindet::Float64
end
SpreadHTCM( winlen, s ) = SpreadHTCM( winlen, s, 1e-12 )
struct DirectionHTCM{T,S} <: HeadTailClassificationMethod
    winlen::Int
    t::T
    s::S
end

headtail_stats( m::SpeedHTCM, e1, e2, center ) = speed_stats( e1, e2, center, m.t, m.s )
headtail_stats( m::SpreadHTCM, e1, e2, center ) = spread_stats( e1, e2, center, m.winlen, m.s, m.mindet )
headtail_stats( m::DirectionHTCM, e1, e2, center ) =  direction_stats(e1, e2, center, m.t, m.winlen, m.s)
