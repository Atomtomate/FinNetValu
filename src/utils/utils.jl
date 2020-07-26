"""
    constantly(val)

Create constant function that returns `val` when called.
"""
constantly(val) = (x...) -> val

"""
    calm(f, n)

Create a new function that behaves like `f()` but returns the same
value for `n` successive calls, i.e. `f` is invoked once for every `n`
calls of `calm(f, n)`.  
"""
function calm(f, n)
    i = 0
    val = f()
    function calmed()
        if i == n
            i = 1
            val = f()
        else
            i += 1
        end
        val
    end
    calmed
end

"""
    bisect(f, l, r, args...; maxiter = 10^2, allow_swap = false)

returns intervall [`l'`,`r'`] with `f(x, args...) == 0` for x ∈ [`l'`,`r'`] and
a function `f` with f > 0 for x' > 0 and f < 0 for x' < 0 in the [`l`, `r`],
using the bisection algorithm. If `allow_swap` is `true`, the order 
of `l` and `r` may be reversed.
The maximum number of bisections can be specified with `maxiter`.
The minimum size of the interval [`l'`,`r'`] can be specified with `isize`.
"""
function bisect(f, l, r, args...; maxiter = 10^12, isize = 10^-12, allow_swap = false)
    if allow_swap
        l, r = (l > r) ? (r, l) : (l, r)
    end
    for i in 1:maxiter
        mid = (r+l)/2
        res = f(mid, args...)
        l = res < 0 ? mid : l
        r = res > 0 ? mid : r
        res == 0 && return [l,r]
        (r - l) < isize && return [l,r]
    end
    error("Maximum number of iterations exceeded")
end
