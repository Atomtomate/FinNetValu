import Base.getindex
import Base.length
import Base.size
import Base.lastindex

"""
    FixedPoint(x)

FixedPoint for Financial Model. Constant wrapper around `Array`,
will potentially hold additional information in future.
"""
struct FixedPoint <: AbstractArray{Float64,1}
    x::Array{Float64,1}
end
const StateVector = Union{FixedPoint, Array}

getindex(fp::FixedPoint, ind...) = fp.x[ind...];
getindex(fp::FixedPoint) = fp.x;
size(fp::FixedPoint) = size(fp.x)
length(fp::FixedPoint) = length(fp.x)
lastindex(fp::FixedPoint) = lastindex(fp.x)

"""
    fixedpoint_naive(mapping, init; maxit=10e8, debug=false)

Find fixedpoint for mapping by simple jocobi iteration.
This is a naive default implementation.
"""
function fixedpoint_naive(mapping, init; maxit=10e8, debug=false)
    xold = init
    xnew = mapping(init)
    it = 1

    while !all(xold .≈ xnew)
        debug && println(xold)
        xold = xnew
        xnew = mapping(xnew)
        it > maxit && throw("Unable to converge fixed point!")
    end
    return FixedPoint(xnew)
end
