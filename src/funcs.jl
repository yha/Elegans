circres(v) = mean(exp.(v.*im))
circmean(v) = angle(circres(v))
circvar(v) = 1 - norm(circres(v))

fix_angle(x) = mod(x+π,2π)-π
m2n(a) = coalesce.(a,NaN)

## Missing Missing methods

Base.atan(::Union{Real,Missing}, ::Union{Real,Missing}) = missing
Base.hypot(::Union{Number,Missing}...) = missing
Base.angle(::Missing) = missing

if VERSION < v"1.2"
    Base.round(::Missing; digits=0, base=0) = missing
    Base.zero(::Type{Missing}) = missing
    Base.oneunit(::Type{Missing}) = missing
end

import StatsBase
(::StatsBase.ECDF)(::Missing) = missing
