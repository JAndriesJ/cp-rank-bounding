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
                H[[i1,i2,j1,j2]]  = i1j1 + i2j2
            end
        end
    end
end

##Mixed examples

# Randomly generated matrices:
""" Generates a sep-state ρ with size: d₁×d₂ and cp-rank: r."""
function genSEPmatrix(d₁,d₂,r)
    A = zeros(Float64, (d₁*d₂,d₁*d₂))
    factor_dict = Dict()
    for ℓ in 1:r
        a_ℓ = rand(Float64, (d₁, 1)) + rand(Float64, (d₁, 1))*im
        b_ℓ = rand(Float64, (d₂, 1)) + rand(Float64, (d₂, 1))*im
        #factor_dict[ℓ]
        A = A + kron(a_ℓ* transpose(conj(a_ℓ)), b_ℓ* transpose(conj(b_ℓ)))
    end
    A = makediagone(A) #is this still relevant?
end


# Utility functions:
