# This script is a colllection of completely positive matrices used for the testing of cp-rank computation.

using LinearAlgebra # For the matrix manipulations and preprocessing.
using Random # for the generation of a random cp Matrix
using Test # Utility module for testing if code works as intended.
Random.seed!(1234)

"""Vector input, output a Circulant matrix of said vector"""
function circmatrix(v)
    hcat([circshift(v, k) for k = -length(v)+1:0]...)
end

"""Normalization function"""
makediagone(A) = diagm(1 ./ sqrt.(diag(A))) * A * diagm( 1 ./ sqrt.(diag(A)))

# Matrices from:  I.M. Bomze, W. Schachinger, R. Ullrich. From seven to eleven: Completely positive matrices with high cp-rank. Linear Algebra and its Applications 459 (2014), 208 – 221.
M7  = circmatrix([531.0, 81, 150, 150, 81, 531, 926]/926)

M7tilde = circmatrix([108, 27, 4, 4, 27, 108, 163.0]/163)

M8tilde = makediagone([541.0 880 363 24 55 11 24 0;
                       880 2007 1496 363 48 22 22 24;
                       363 1496 2223 1452 363 24 22 11;
                       24 363 1452 2325 1584 363 48 55;
                       55 48 363 1584 2325 1452 363 24;
                       11 22 24 363 1452 2223 1496 363;
                       24 22 22 48 363 1496 2007 880;
                       0 24 11 55 24 363 880 541])

M9tilde = makediagone([2548 1628 363 60 55 55 60 363 1628;
                       1628 2548 1628 363 60 55 55 60 363;
                       363 1628 2483 1562 363 42 22 55 60;
                       60 363 1562 2476 1628 363 42 55 55;
                       55 60 363 1628 2548 1628 363 60 55;
                       55 55 42 363 1628 2476 1562 363 60;
                       60 55 22 42 363 1562 2483 1628 363;
                       363 60 55 55 60 363 1628 2548 1628;
                       1628 363 60 55 55 60 363 1628 2548])

M11tilde =      makediagone([781 0 72 36 228 320 240 228 36 96 0;
                             0 845 0 96 36 228 320 320 228 36 96;
                            72 0 827 0 72 36 198 320 320 198 36;
                            36 96 0 845 0 96 36 228 320 320 228;
                            228 36 72 0 781 0 96 36 228 240 320;
                            320 228 36 96 0 845 0 96 36 228 320;
                            240 320 198 36 96 0 745 0 96 36 228;
                            228 320 320 228 36 96 0 845 0 96 36;
                            36 228 320 320 228 36 96 0 845 0 96;
                            96 36 198 320 240 228 36 96 0 745 0;
                            0 96 36 228 320 320 228 36 96 0 845])


# Randomly generated matrices:
""" Generates a cp-matrix with size: n and cp-rank at most r."""
function genCPmatrix(n,r)
    A = zeros(Float64, (n,n))
    factor_dict = Dict()
    for ℓ in 1:r
        a_ℓ = rand(Float64, (n, 1))
        #factor_dict[ℓ]
        A = A + a_ℓ* transpose(a_ℓ)
    end
    A = makediagone(A)
end
@testset "genCPmatrix" begin
A = genCPmatrix(4,2)
@test all(A .>= 0)
eigs = eigvals(A)
@test all(eigs .+ 1.0e-14  .>= 0)
end
