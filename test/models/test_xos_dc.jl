using ForwardDiff
using LinearAlgebra
using SparseArrays


@testset "xos default cost general" begin
    N  = 4
    Mˢ = [0.0 0.2 0.3 0.1; 0.2 0.0 0.2 0.1; 0.1 0.1 0.0 0.3; 0.1 0.1 0.1 0.0]
    Mᵈ = [0.0 0.0 0.1 0.0; 0.0 0.0 0.0 0.1; 0.1 0.0 0.0 0.1; 0.0 0.1 0.0 0.0]
    Mᵉ = I
    d1  = [0.0, 0.0, 0.0, 0.0]
    d2  = [0.8, 0.8, 0.8, 0.8]
    α1 = I
    α2 = [0.5, 0.5, 0.5, 0.5]
    net1 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α1)
    net2 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α2)
    net3 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α1)
    net4 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α2)
    a = [2.0, 0.5, 0.6, 0.6]
    x1 = FinNetValu.fixvalue(net1, a)
    println("Solvent: ", round.(FinNetValu.solvent(net1, x1), digits=3))
    println("Fixed Point: ", round.(x1, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net1, x1, a), digits=3))
    x2 = FinNetValu.fixvalue(net2, a)
    println("Solvent: ", round.(FinNetValu.solvent(net2, x2), digits=3))
    println("Fixed Point: ", round.(x2, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net2, x2, a), digits=3))
    x3 = FinNetValu.fixvalue(net3, a)
    println("Solvent: ", round.(FinNetValu.solvent(net3, x3), digits=3))
    println("Fixed Point: ", round.(x3, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net3, x3, a), digits=3))
    x4 = FinNetValu.fixvalue(net4, a)
    println("Solvent: ", round.(FinNetValu.solvent(net4, x4), digits=3))
    println("Fixed Point: ", round.(x4, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net4, x4, a), digits=3))
    # TODO: auf papier nachrechnen und pruefen
end
