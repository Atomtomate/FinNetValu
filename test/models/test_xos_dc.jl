using ForwardDiff
using LinearAlgebra
using SparseArrays

@testset "xos default" begin
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
    net2 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α0)
    x2 = FinNetValu.fixvalue(net2, a)
    net3 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α0)
    x3 = FinNetValu.fixvalue(net3, a)

    net4 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d0, α1)
    x4 = FinNetValu.fixvalue(net4, a)
    net5 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d1, α1)
    x5 = FinNetValu.fixvalue(net5, a)
    net6 = FinNetValu.XOSModel(Mˢ, Mᵈ, Mᵉ, d2, α1)
    x6 = FinNetValu.fixvalue(net6, a)

@testset "xos default fixpoint" begin
    @test all(x1 .≈ [16/11, 15/11, 0, 0]) # (1 - Mˢ)^(-1) * (1,1)
    @test all(x2 .≈ [0, 0, 5/3, 4/3])     # (1 - Mᵈ)^(-1) * (1,1)
    @test all(x3 .≈ [12/7, 0, 0, 10/7])     # (1 - Mᵈ)^(-1) * (1,1)

    @test all(x4 .≈ [16/11, 15/11, 0, 0]) # (1 - Mˢ)^(-1) * (1,1)
    @test all(x5 .≈ 1/10 .* [0, 0, 5/3, 4/3])     # (1 - Mᵈ)^(-1) * (1/10,1/10)
    @test all(x6 .≈ [6/5, 0, 0, 2/5])     # (1 - Mᵈ)^(-1) * (1,1)
end

@testset "xos multiple fixed points" begin
    Mᵉ₂ = I
    Mˢ₂ = zeros(2,2)
    α0₂ = [1/2, 1/2]
    L =  [0 2.2; 1/(1-α0₂[1]) 0]
    Lbar = sum(L, dims=[2])[:,1]
    Π₂ = L ./ Lbar
    Mᵈ₂ = transpose(Π₂)
    a₂  = [1.0, 1.0]
    d0₂ = Lbar
    net1₂ = FinNetValu.XOSModel(Mˢ₂, Mᵈ₂, Mᵉ₂, d0₂, α0₂)
end

@testset "xos default N = 4" begin
    Mˢ₃ = [0.0 0.2 0.3 0.1; 0.2 0.0 0.2 0.1; 0.1 0.1 0.0 0.3; 0.1 0.1 0.1 0.0]
    Mᵈ₃ = [0.0 0.0 0.1 0.0; 0.0 0.0 0.0 0.1; 0.1 0.0 0.0 0.1; 0.0 0.1 0.0 0.0]
    Mᵉ₃ = I
    d1₃  = [0.0, 0.0, 0.0, 0.0]
    d2₃  = [0.8, 0.8, 0.8, 0.8]
    α1₃ = I
    α2₃ = [0.5, 0.5, 0.5, 0.5]
    net1₃ = FinNetValu.XOSModel(Mˢ₃, Mᵈ₃, Mᵉ₃, d1₃, α1₃)
    net2₃ = FinNetValu.XOSModel(Mˢ₃, Mᵈ₃, Mᵉ₃, d1₃, α2₃)
    net3₃ = FinNetValu.XOSModel(Mˢ₃, Mᵈ₃, Mᵉ₃, d2₃, α1₃)
    net4₃ = FinNetValu.XOSModel(Mˢ₃, Mᵈ₃, Mᵉ₃, d2₃, α2₃)
    a₃ = [2.0, 0.5, 0.6, 0.6]
    x1₃ = FinNetValu.fixvalue(net1₃, a₃)
    println("Solvent: ", round.(FinNetValu.solvent(net1₃, x1₃), digits=3))
    println("Fixed Point: ", round.(x1₃, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net1₃, x1₃, a₃), digits=3))
    x2₃ = FinNetValu.fixvalue(net2₃, a₃)
    println("Solvent: ", round.(FinNetValu.solvent(net2₃, x2₃), digits=3))
    println("Fixed Point: ", round.(x2₃, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net2₃, x2₃, a₃), digits=3))
    x3₃ = FinNetValu.fixvalue(net3₃, a₃)
    println("Solvent: ", round.(FinNetValu.solvent(net3₃, x3₃), digits=3))
    println("Fixed Point: ", round.(x3₃, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net3₃, x3₃, a₃), digits=3))
    x4₃ = FinNetValu.fixvalue(net4₃, a₃)
    println("Solvent: ", round.(FinNetValu.solvent(net4₃, x4₃), digits=3))
    println("Fixed Point: ", round.(x4₃, digits=3))
    println("Valuation: ", round.(FinNetValu.valuation(net4₃, x4₃, a₃), digits=3))
    # TODO: auf papier nachrechnen und pruefen
end

@testset "default boundaries" begin
    #TODO: calculate positions of boundaries along one dimension (bisection)
    #TODO: calculate all 4 bondaries for N=2
    #TODO: recursive algorithm for N > 2
end

@testset "default gaps" begin
    #TODO: calculate loss for all default boundaries
end
end
