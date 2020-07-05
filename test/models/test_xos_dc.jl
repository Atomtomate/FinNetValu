using ForwardDiff
using LinearAlgebra
using SparseArrays

@testset "xos default fixpoint" begin
    Mˢ = [0.0 1/3; 1/4 0.0]
    Mᵈ = [0.0 1/2; 1/5 0.0]
    Mᵉ = I
    a  = [1.0, 1.0]
    d0 = [0, 0]
    d1 = [5, 5]
    d2 = [0, 5]
    α0 = [1, 1]
    α1 = [1/10, 1/10]

    net1 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d0, α0)
    x1 = FinNetValu.fixvalue(net1, a)
    @test all(x1 .≈ [16/11, 15/11, 0, 0]) # (1 - Mˢ)^(-1) * (1,1)
    net2 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α0)
    x2 = FinNetValu.fixvalue(net2, a)
    @test all(x2 .≈ [0, 0, 5/3, 4/3])     # (1 - Mᵈ)^(-1) * (1,1)
    net3 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α0)
    x3 = FinNetValu.fixvalue(net3, a)
    @test all(x3 .≈ [12/7, 0, 0, 10/7])     # (1 - Mᵈ)^(-1) * (1,1)

    net4 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d0, α1)
    x4 = FinNetValu.fixvalue(net4, a)
    @test all(x4 .≈ [16/11, 15/11, 0, 0]) # (1 - Mˢ)^(-1) * (1,1)
    net5 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α1)
    x5 = FinNetValu.fixvalue(net5, a)
    @test all(x5 .≈ (1/10) .* [0, 0, 5/3, 4/3])     # (1 - Mᵈ)^(-1) * (1/10,1/10)
    net6 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α1)
    x6 = FinNetValu.fixvalue(net6, a)
    @test all(x6 .≈ [6/5, 0, 0, 2/5])     # (1 - Mᵈ)^(-1) * (1,1)

    println(x6)
end

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
