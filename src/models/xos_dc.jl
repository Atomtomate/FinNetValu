
function loss_db(net, aₗ, aᵣ) 
    xₗ = FinNetValu.fixvalue(net, aₗ)
    xᵣ = FinNetValu.fixvalue(net, aᵣ)
    [sum(abs.(FinNetValu.equity(xₗ) .- FinNetValu.equity(xᵣ))), 
     sum(abs.(FinNetValu.debt(xₗ) .- FinNetValu.debt(xᵣ)))]
end

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

function find_db(net, aₗ, aᵣ, ζ, b; atol = 10^-10)
    aint = copy(aₗ)
    function fi(ai)
        aint[b] = ai
        res = 2*sum(FinNetValu.solvent(net, FinNetValu.fixvalue(net, aint)) .- ζ) .- 1
        return res
    end
    return bisect(fi, aₗ[b], aᵣ[b])
end

function map_db(f, net, a)
    x = FinNetValu.fixvalue(net, a)
    defaulting_banks = findall(xi -> xi == 0,FinNetValu.solvent(net, x))
    losses = zeros(length(defaulting_banks),2)
    ai = copy(a)

        
    for (bi,b) in enumerate(defaulting_banks)
        ai_l = copy(ai)
        ai_r = copy(ai)
        ζ_l = FinNetValu.solvent(net, FinNetValu.fixvalue(net, ai))
        a_r = find_r(net, ai, ζ_l, b)
        ai_l[b], ai_r[b] = find_db(net, ai, a_r, ζ_l, b)
        losses[b, :] = f(net, ai_l, ai_r)
        ai = a_r
    end
    return losses
end
