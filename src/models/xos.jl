import Base.show
using LinearAlgebra

const MatrixType{T} = Union{AbstractMatrix{T1}, UniformScaling{T2}} where {T1<:Real, T2<:Real}
const VectorType{T} = Union{AbstractVector{T1}, UniformScaling{T2}} where {T1<:Real, T2<:Real}

"""
    XOSModel(N, Mˢ, Mᵈ, Mᵉ, d)

Financial network of `N` firms with investment portfolios `Mˢ`, `Mᵈ`
and `Mᵉ` of holding fractions in counterparties equity, debt and
external assets respectively. Nominal debt `d` is due at maturity.
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

function valuation!(y, net::XOSModel, x, a, dc=nothing)
    dc = dc === nothing ? defaultcosts(net, x) : dc
    tmp =  net.Mˢ * equityview(net, x) .+ net.Mᵈ * debtview(net, x)
    equityview(net, y) .= max.(zero(eltype(x)), net.Mᵉ * a .+ tmp .- nominaldebt(net))
    debtview(net, y)   .= min.(nominaldebt(net), dc .* (net.Mᵉ * a) .+ tmp)
end

function valuation(net::XOSModel, x, a, dc=nothing)
    dc = dc === nothing ? defaultcosts(net, x) : dc
    tmp = net.Mˢ * equityview(net, x) .+ net.Mᵈ * debtview(net, x)
    vcat(max.(zero(eltype(x)), (net.Mᵉ * a) .+ tmp .- nominaldebt(net)),
         min.(nominaldebt(net), dc .* (net.Mᵉ * a) .+ tmp))
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

function solvent(net::XOSModel, x)
    equityview(net, x) .> zero(eltype(x))
end

init(net::XOSModel) = repeat(numfirms(net) .* nominaldebt(net), 2)

##########################
# Model specific methods #
##########################
function fixedpoint_naive(mapping, init, maxit=10e8)
    xold = init
    xnew = mapping(init)
    it = 1

    while !all(xold .≈ xnew)
        println(xold)
        xold = xnew
        xnew = mapping(xnew)
        it > maxit && throw("Unable to converge fixed point!")
    end
    return xnew
end

"""
    TODO: this is a primitive test in order to check for convergence of
    of the double fixed point iteration to an arbitrary default cost vector.
    This needs to be extended to allow for default cost functions.
"""
function fixvalue(net::XOSModel, a; dc_it = 5 , m = 0, kwargs...)
    x_i = init(net)
    ones = one(eltype(defaultcosts(net, x_i)))
    #TODO: in general, this should be a function, giving the next dc_i, converged is it does not change
    x_0 = fixedpoint_naive(x-> valuation(net, x, a, fill(1.0,   length(a))), x_i)
    x   = fixedpoint_naive(x-> valuation(net, x, a, net.αᵉ), x_0)
    #= for dc_i in (ones .- i*(ones .- net.αᵉ)/dc_it for i in 0:dc_it) =#
    #=     x_i = fixedpoint_naive( x-> valuation(net, x, a, dc_i), x_i) =#
    #=     println("it: ", round.(x_i,digits=3), " for dc = ", round.(dc_i, digits=3)) =#
    #= end =#
    return x
end


"""
    equityview(net, x)

View the equity part of `x` which can be a state vector of Jacobian
matrix.
"""
function equityview end

"""
    debtview(net, x)

View the debt part of `x` which can be a state vector of Jacobian
matrix.
"""
function debtview end

equityview(net::XOSModel, x::AbstractVector) = view(x, 1:numfirms(net))
#equityview(net::XOSModel, x::AbstractMatrix) = view(x, 1:numfirms(net), :)

debtview(net::XOSModel, x::AbstractVector) = begin N = numfirms(net); view(x, (N+1):(2*N)) end
debtview(net::XOSModel, x::AbstractMatrix) = begin N = numfirms(net); view(x, (N+1):(2*N), :) end

@inline defaultcosts(net::XOSModel, x, dc=net.αᵉ) = dc .* ( .! solvent(net, x))
@inline nominaldebt(net::XOSModel) = net.d
