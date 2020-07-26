"""
    find_r(net, a, ζ, b; maxiter = 10^5)

Returns `a'` for wich firm `b` becomes solvent, assuming
`b` is not solvent under the fixed point solution for `net` and `a`.
`ζ` is the solvency vector for the initial fixed point.
Note, that this is NOT the smallest possible `a'`, but it is guaranteed,
that no other firm changes its solvency status.
"""
function find_r(net, a, ζ, b; maxiter = 10^5)
    aint = copy(a)
    
    function f(ai, a)
        aint[b] = ai
        sum(ζ .- FinNetValu.solvent(net, FinNetValu.fixvalue(net, aint))) + 1
    end    
    l,r = bisect(f, a[b], net.N * net.d[b], a)
    aint[b] = (r+l)/2
    return aint
end

"""
    find_db(net, aₗ, aᵣ, ζ, b; atol = 10^-10)

Returns an interval for the default boundary of firm `b`.
[`aₗ`, `aᵣ`] specifies the initial search intervall in which only
firm `b` changes solvency status. This interval can be obtained
using `find_r`.
`atol` specifies the maximum width of the interval.
`ζ` is the solvency vector for the initial fixed point.
"""
function find_db(net, aₗ, aᵣ, ζ, b; atol = 10^-10)
    aint = copy(aₗ)
    function fi(ai)
        aint[b] = ai
        res = 2*sum(FinNetValu.solvent(net, FinNetValu.fixvalue(net, aint)) .- ζ) .- 1
        return res
    end
    return bisect(fi, aₗ[b], aᵣ[b])
end

"""
    map_db(f, net, a)

Maps a function f(net::FinancialModel, a1, a2) over all default boundary 
crossings between the fixed point corresponding to `net` and `a` up to 
`a'` for which all firms are solvent. 
"""
function map_db(f, net, a)
    x = fixvalue(net, a)
    defaulting_banks = findall(xi -> xi == 0, solvent(net, x))
    losses = zeros(length(defaulting_banks),2)
    ai = copy(a)
        
    for (bi,b) in enumerate(defaulting_banks)
        ai_l = copy(ai)
        ai_r = copy(ai)
        ζ_l = FinNetValu.solvent(net, FinNetValu.fixvalue(net, ai))
        a_r = find_r(net, ai, ζ_l, b)
        ai_l[b], ai_r[b] = find_db(net, ai, a_r, ζ_l, b)
        losses[bi, :] = f(net, ai_l, ai_r)
        ai = a_r
    end
    return losses
end
