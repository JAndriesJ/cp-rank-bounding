using LinearAlgebra # For the matrix manipulations and preprocessing.
using Random # for the generation of a random cp Matrix
using Test # Utility module for testing if code works as intended.

using Kronecker

Random.seed!(1234)

## Pure examples

#Example 3.7. Consider the Hankel tensor H ∈ C^[2,2] in [NieYe19] such that
H  = Dict()
n = 2
for i1 in 1:n
    for i2 in 1:n
        for j1 in 1:n
            for j2 in 1:n
                H[[i1,i2,j1,j2]]  = i1 + i2 + j1 + j2
            end
        end
    end
end

#Example 3.7. Consider the Hankel tensor H ∈ C^[2,2] in [NieYe19] such that
H  = Dict()
n = 3
for i1 in 1:n
    for i2 in 1:n
        for j1 in 1:n
            for j2 in 1:n
                H[[i1,i2,j1,j2]]  = i1*j1 + i2*j2
            end
        end
    end
end

##Mixed examples


# Randomly generated matrices:
""" Generates a sep-matrix with size: d₁*d₂ and Sep-rank at most r."""
function genRandGSEPmatrix(d₁,d₂,r,isComp = false)
    ρ = zeros(Float64, (d₁*d₂,d₁*d₂))
    # factor_dict = Dict()
    a_ℓ = Dict()
    b_ℓ = Dict()
    for ℓ in 1:r
        if isComp
            a_ℓtemp = randn(Float64, (d₁, 1)) + randn(Float64, (d₁, 1))*im
            b_ℓtemp = randn(Float64, (d₂, 1)) + randn(Float64, (d₂, 1))*im
        else
            a_ℓtemp = randn(Float64, (d₁, 1))
            b_ℓtemp = randn(Float64, (d₂, 1))
        end
        a_ℓ[ℓ] = a_ℓtemp / norm(a_ℓtemp)
        b_ℓ[ℓ] = b_ℓtemp / norm(b_ℓtemp)
        # factor_dict[ℓ]
        ρ = ρ + kron( a_ℓ[ℓ]* transpose(a_ℓ[ℓ]),  b_ℓ[ℓ]* transpose(b_ℓ[ℓ]))
    end
    return ρ
end
@testset "genRandSEPmatrix" begin
ρ = genRandGSEPmatrix(3,3,3)
eigs = eigvals(ρ)
@test all(eigs .+ 1.0e-14  .>= 0)
end


# Utility functions:
