import Base.show
using LinearAlgebra

const MatrixType{T} = Union{AbstractMatrix{T1}, UniformScaling{T2}} where {T1<:Real, T2<:Real}
const VectorType{T} = Union{AbstractVector{T1}, UniformScaling{T2}} where {T1<:Real, T2<:Real}

"""
    XOSModel(Mˢ, Mᵈ, Mᵉ, d, αᵉ)

Financial network of `N` firms with investment portfolios `Mˢ`, `Mᵈ`
and `Mᵉ` of holding fractions in counterparties equity, debt and
external assets respectively. Nominal debt `d` is due at maturity.
The default costs are given by `αᵉ`.
"""
struct XOSModel <: FinancialModel
    N::Int64
    Mˢ::MatrixType
    Mᵈ::MatrixType
    Mᵉ::MatrixType
    d::VectorType
    αᵉ::AbstractVector

    function XOSModel(Mˢ::MatrixType, Mᵈ::MatrixType, Mᵉ::MatrixType, d::VectorType, αᵉ::VectorType)
        @assert isleft_substochastic(Mˢ)
        @assert isleft_substochastic(Mᵈ)
        @assert all(d .>= 0)
        ld = length(collect(d))
        if typeof(Mˢ) <: Array
            @assert ld == size(Mˢ,1)
            @assert ld == size(Mˢ,2)
        end
        if typeof(Mᵈ) <: Array
            @assert ld == size(Mᵈ,1)
            @assert ld == size(Mᵈ,2)
        end
        if typeof(Mᵉ) <: Array
            @assert ld == size(Mᵉ,1)
            @assert ld == size(Mᵉ,2)
        end
        αᵉ_internal = αᵉ isa UniformScaling ? fill(αᵉ.λ,ld) : αᵉ
        new(ld, Mˢ, Mᵈ, Mᵉ, d, αᵉ_internal)
    end

    function XOSModel(Mˢ::MatrixType{T1,T2}, Mᵈ::MatrixType{T1,T2}, Mᵉ::MatrixType{T1,T2},
                      d::AbstractVector{T1}) where {T1<:Real, T2<:Real}
        XOSModel(Mˢ, Mᵈ, Mᵉ, d, I)               
    end
end

function show(io::IO, net::XOSModel)
    if any(net.Mˢ .> 0) & any(net.Mᵈ .> 0)
        msg = "equity and debt"
    elseif any(net.Mˢ .> 0)
        msg = "equity"
    elseif any(net.Mᵈ .> 0)
        msg = "debt"
    else
        msg = "no"
    end
    print(io, "XOS model of N = ",
          numfirms(net), " firms with ",
          msg, " cross holdings.")
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::XOSModel) = net.N

function valuation(net::XOSModel, x::StateVector, a, dc=nothing)
    tmp = net.Mˢ * equity(x) .+ net.Mᵈ * debt(x)
    vcat(max.(zero(eltype(x)), (net.Mᵉ * a) .+ tmp .- nominaldebt(net)),
         min.(nominaldebt(net), defaultcosts(net, x, dc) .* (net.Mᵉ * a) .+ tmp))
end

function fixjacobian(net::XOSModel, a, x = fixvalue(net, a))
    ## Uses analytical formulas for speedup
    ξ = solvent(net, x)
    eins = one(eltype(ξ))
    dVdx = vcat(hcat(Diagonal(ξ) * net.Mˢ, Diagonal(ξ) * net.Mᵈ),
                hcat(Diagonal(eins .- ξ) * net.Mˢ, Diagonal(eins .- ξ) * net.Mᵈ))
    dVda = vcat(Diagonal(ξ), Diagonal(1.0 .- ξ)) * net.Mᵉ
    (I - dVdx) \ Matrix(dVda) ## Note: RHS needs to be dense
end

function solvent(net::XOSModel, x::StateVector)
    equity(x) .> zero(eltype(x))
end

init(net::XOSModel) = repeat(numfirms(net) .* nominaldebt(net), 2)

##########################
# IntervalBox Functions  #
##########################
function valuation(net::XOSModel, x::IntervalBox, a, dc=nothing)
    N = numfirms(net)
    tmp = net.Mˢ * x[1:N] .+ net.Mᵈ * x[N+1:2*N]
    return SVector(max.(zero(eltype(x)), (net.Mᵉ * a) .+ tmp .- net.d)...,
                   min.(net.d, defaultcosts(net, x, dc).* (net.Mᵉ * a) .+ tmp)...)
end

function solvent(net::XOSModel, x::IntervalBox)
    x[1:numfirms(net)] .> zero(eltype(x))
end



##########################
# Model specific methods #
##########################

"""
    TODO: this is a primitive test in order to check for convergence of
    of the double fixed point iteration to an arbitrary default cost vector.
    This needs to be extended to allow for default cost functions.
"""
function fixvalue(net::XOSModel, a; dc_it = 5 , m = 0, kwargs...)
    x_i = init(net)
    ones = one(eltype(defaultcosts(net, x_i)))
    #TODO: in general, this should be a function, giving the next dc_i, converged is it does not change
    x_0 = fixedpoint_naive(x-> valuation(net, x, a, 1.0), x_i, debug=false)
    x   = fixedpoint_naive(x-> valuation(net, x, a, net.αᵉ), x_0, debug=false)
    #= for dc_i in (ones .- i*(ones .- net.αᵉ)/dc_it for i in 0:dc_it) =#
    #=     x_i = fixedpoint_naive( x-> valuation(net, x, a, dc_i), x_i) =#
    #=     println("it: ", round.(x_i,digits=3), " for dc = ", round.(dc_i, digits=3)) =#
    #= end =#
    return x
end

equity(x::StateVector) = begin N = floor(Int64,length(x)/2); x[1:N] end
debt(x::StateVector) = begin N = floor(Int64,length(x)/2); x[N+1:end] end
equityview(net::XOSModel, x::Array) = view(x, 1:numfirms(net))
debtview(net::XOSModel, x::Array) = begin N = numfirms(net); view(x, (N+1):(2*N)) end


function defaultcosts(net::XOSModel, x, dc=nothing)
    ζ = solvent(net, x)
    dc = dc === nothing ? net.αᵉ : dc 
    return dc .* ( .!ζ).+ 1.0 .* ζ
end
nominaldebt(net::XOSModel) = net.d


############################
# Default Boundary Related #
############################
function loss_db(net, aₗ, aᵣ) 
    xₗ = FinNetValu.fixvalue(net, aₗ)
    xᵣ = FinNetValu.fixvalue(net, aᵣ)
    [sum(abs.(FinNetValu.equity(xₗ) .- FinNetValu.equity(xᵣ))), 
     sum(abs.(FinNetValu.debt(xₗ) .- FinNetValu.debt(xᵣ)))]
end
