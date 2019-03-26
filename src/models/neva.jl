import Base.show

"""
    NEVAModel(Lᵉ, L, 𝕍ᵉ, 𝕍)

Financial network model with nominal external `Lᵉ` and internal
liabilities `L`. Nominal values of external assets and equity are
weighted by the valuation functions `𝕍ᵉ` and `𝕍` respectively.
"""
struct NEVAModel <: FinancialModel
    name
    N
    A
    l
    𝕍ᵉ
    𝕍

    """
        NEVAModel(N, Lᵉ, L, 𝕍ᵉ, 𝕍)

    Construct NEVA model with `N` firms, external `Lᵉ` and internal
    liabilities `L`. Values of external assets and liabilities of
    counterparties are adjusted by the valuation functions `𝕍ᵉ` and
    `𝕍` respectively."""

    function NEVAModel(name::String, N::Integer, Lᵉ, L, 𝕍ᵉ::Function, 𝕍::Function)
        @assert all(Lᵉ .>= 0)
        @assert all(L .>= 0)
        A = copy(L')
        l = rowsums(L) .+ Lᵉ
        new(name, N, A, l, 𝕍ᵉ, 𝕍)
    end
end

function show(io::IO, net::NEVAModel)
    print(io, net.name, " model with N = ", numfirms(net), " firms.")
end

############################
# Convenience constructors #
############################

function NEVAModel(name::String, Lᵉ::AbstractVector, L, 𝕍ᵉ::Function, 𝕍::Function)
    NEVAModel(name, length(Lᵉ), Lᵉ, L, 𝕍ᵉ, 𝕍)
end

function NEVAModel(name::String, Lᵉ, L::AbstractMatrix, 𝕍ᵉ::Function, 𝕍::Function)
    @assert size(L, 1) == size(L, 2)
    NEVAModel(name, size(L, 1), Lᵉ, L, 𝕍ᵉ, 𝕍)
end

function NEVAModel(name::String, Lᵉ::AbstractVector, L::AbstractMatrix, 𝕍ᵉ::Function, 𝕍::Function)
    @assert length(Lᵉ) == size(L, 1) == size(L, 2)
    NEVAModel(name, length(Lᵉ), Lᵉ, L, 𝕍ᵉ, 𝕍)
end

##############################################
# Implementation of FinancialModel interface #
##############################################

numfirms(net::NEVAModel) = net.N

function valuation!(y, net::NEVAModel, x, a)
    y .= a .* net.𝕍ᵉ(net, x, a) .+ rowsums(net.A .* net.𝕍(net, x, a)) .- net.l
end

function valuation(net::NEVAModel, x, a)
    a .* net.𝕍ᵉ(net, x, a) .+ rowsums(net.A .* net.𝕍(net, x, a)) .- net.l
end

function solvent(net::NEVAModel, x)
    x .> zero(eltype(x))
end

function init(net::NEVAModel, a)
    ## Initialize between boundaries m <= M
    a .- net.l
end

##########################
# Model specific methods #
##########################

bookequity(net::NEVAModel, a) = a .+ rowsums(net.A) .- net.l

##########################################
# Constructors for different models from #
# arxiv???                               #
##########################################

"""
    EisenbergNoeModel(Lᵉ, L)

Creates an instance of the NEVAModel with valuation functions

```math
\\begin{align}
\\mathbb{V}^e_i(E_i) &= 1 \\quad \\forall i \\\\
\\mathbb{V}_{ij}(E_j) &= \\unicode{x1D7D9}_{E_j \\geq 0} + \\left(\\frac{E_j + \\bar{p}_j}{\\bar{p}_j}\\right)^+ \\unicode{x1D7D9}_{E_j < 0} \\quad \\forall i, j
\\end{align}
```

where ``\\bar{p}_j = \\sum_k L_{jk}``.

This valuation was shown to correspond to the model by Eisenberg &
Noe.
"""
function EisenbergNoeModel(Lᵉ::AbstractVector, L::AbstractMatrix)
    pbar = vec(sum(L; dims = 2))
    function val(net, e, a)
        funEN(e, pbar) = if (e > 0) 1. else max((e + pbar) / pbar, 0.) end
        # Note: rowvector gets broadcasted correctly as 𝕍(Eⱼ)
        transpose(funEN.(e, pbar))
    end
    NEVAModel("Eisenberg & Noe",
              Lᵉ,
              L,
              constantly(one(eltype(L))),
              val)
end

function FurfineModel(Lᵉ::AbstractVector, L::AbstractMatrix, R::Real)
    @assert 0 <= R <= 1
    function val(net, e, a)
        funF(e) = if (e > 0) 1. else R end
        transpose(funF.(e))
    end
    NEVAModel("Furfine",
              Lᵉ,
              L,
              constantly(one(eltype(L))),
              val)
end

function LinearDebtRankModel(Lᵉ::AbstractVector, L::AbstractMatrix)
    function val(net, e, a)
        funLR(e, ebook) = if (e > ebook) 1. elseif (e > 0) e / ebook else 0. end
        transpose(funLR.(e, bookequity(net, a)))
    end
    NEVAModel("Linear Debt Rank",
              Lᵉ,
              L,
              constantly(one(eltype(L))),
              val)
end
