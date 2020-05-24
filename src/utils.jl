## group_pairwise

import IterTools

struct GroupPairwise{I,F}
    pred::F
    xs::I
end
Base.IteratorSize(::Type{<:GroupPairwise}) = Base.SizeUnknown()
Base.eltype(::Type{<:GroupPairwise{I}}) where {I} = Vector{eltype(I)}

group_pairwise(f,xs) = GroupPairwise(f,xs)

function Base.iterate( it::GroupPairwise{I,F}, state=nothing ) where {I,F}
    if state === nothing
        val, xs_state = IterTools.@ifsomething iterate(it.xs)
    else
        done, val, xs_state = state
        done && return nothing
    end
    values = Vector{eltype(I)}()
    push!( values, val )

    while true
        next = iterate( it.xs, xs_state )
        next === nothing && return values, (true, nothing, nothing)

        prev_val = val
        val, xs_state = next
        if it.pred( prev_val, val )
            push!( values, val )
        else
            break
        end
    end

    return values, (false, val, xs_state)
end

## returning exception

errfilter(f, err_pred, g_ok, g_err) = function(args...)
    errs = Iterators.filter(err_pred, args)
    isempty(errs) ? f((g_ok(a) for a in args)...) : g_err(errs)
end

passex(f) = errfilter(f, x -> x isa Exception, identity, first)
missex(f) = errfilter(f, x -> x isa Exception, identity, _->missing)

try_return(f) = try f(); catch e; e end
trying(f) = (args...) -> try_return(()->f(args...))
# passex(f,x) = x isa Exception ? x : f(x)
# missex(f,x) = x isa Exception ? missing : f(x)
# passex(f) = x -> passex(f,x)
# missex(f) = x -> missex(f,x)

skipex(v) = (x for x in v if !(x isa Exception))

# using ResultTypes
# passerror(f) = errfilter(f, ResultTypes.iserror, ResultTypes.unwrap, first)
# passerror(f) = function (args...)
#     errs = Iterators.filter(ResultTypes.iserror,args)
#     isempty(errs) ? f((unwrap(a) for a in args)...) : first(errs)
# end

## Spreading non-missing values

function spread!(out,in,mask,miss=missing)
    i = firstindex(in)
    for k in eachindex(out)
        if mask[k]
            out[k] = in[i]
            i += 1
        else
            out[k] = miss
        end
    end
    out
end
spread(x,mask,miss=missing) = spread!(similar(mask,Union{eltype(x),typeof(miss)}),x,mask,miss)
#spread_nan(x,mask) = spread(x,mask,NaN)
