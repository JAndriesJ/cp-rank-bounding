#Do you have a moment? (Goal: implement xi_cp of the FOCM paper. Eventually add the G-constraints)
using Test # Utility module for testing if code works as intended.
using LinearAlgebra # For the matrix manipulations and preprocessing.

"""From integers  n and t generate multiindices of N_t^n or N_=t^n  in Lexicographical order."""
function GenMon(n,t, degLEq = true)
    if t == 0
        return zeros(Int32,1,n)
    else
        temp = GenMon(n,t-1,degLEq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if degLEq # For monomials of degree less than t
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end #Has a quadratic loss in runtime because for degree t we compute degree 1 t times...
@testset "GenMon" begin
    MonBase = GenMon(2,3,false)
    @test MonBase ==  [3.0  0.0;
                       2.0  1.0;
                       1.0  2.0;
                       0.0  3.0]

    MonBase = GenMon(2,4,false)
    @test MonBase ==  [4  0;
                       3  1;
                       2  2;
                       1  3;
                       0  4]
end


"""Input is a multi-index array N_t^n and α ∈ N_t^n, output is the row index of the
 monomial with power equal to the multi-index."""
function GetMonIndex(M,α)
    try
        return findall(Bool[transpose(M[i,:]) == α for i =1:size(M,1)])[1] #For some reason the i-th row of an array is returned as a column vector (the transpose remedies this)
    catch
        println("There is no such entry.")
    end
end
@testset "GenMon" begin
    MonBase =  GenMon(2,4,false)
    α = [2  2]
    @test  GetMonIndex(MonBase, α) == 3

    MonBase = GenMon(2,3)
    α = [2.0 1.0]
    @test  GetMonIndex(MonBase, α) == 8
end


"""Convert a sub-array pointer to an array."""
function makeArray(thing)
    α = copy(thing)
    n = length(α)
    α = reshape(α,1,n)
    return α
end


"""Dictionary: keys γ ∈ N_2t^n, values are indeces in Moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ"""
function GenMomMatDict(n,t)
    @assert iseven(t) "$t is an odd number!"
    MonBase = GenMon(n,t/2.0)
    MonMat = Dict()
    for α_slice in eachrow(MonBase)
        α = makeArray(α_slice)
        MonRowInd = GetMonIndex(MonBase, α)
        for β_slice in eachrow(MonBase)
            β = makeArray(β_slice)
            MonColInd = GetMonIndex(MonBase, β)
            γ = α + β
            if haskey(MonMat, γ)
                # MonMat[tempEntry]  = push!(MonMat[γ], [α, β]) #  Index by exponents
                MonMat[γ]  = push!(MonMat[γ], (MonRowInd, MonColInd)) # Index by position in Monmial vector.

            else
                # MonMat[tempEntry]  = [[α, β]] #  Index by exponents
                MonMat[γ]  = [(MonRowInd, MonColInd)]
            end
        end
    end
    return MonMat
end
@testset "GenMomMatDict" begin
    n = 2
    t = 2
    MonExp = GenMomMatDict(n,t*2)
    MonBase = GenMon(n,t)
    for key in keys(MonExp)
        for val in MonExp[key]
            @test   makeArray(key) == makeArray(MonBase[val[1],:] + MonBase[val[2],:])
        end
    end
end

"""Dictionary: keys γ ∈ N_2t^n, values are indeces in Moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ"""
function GenMomMatDictComplex(n,t)
    @assert iseven(t) "$t is an odd number!"
    MonBase = GenMon(n,t/2.0)
    MonMat = Dict()
    for α_slice in eachrow(MonBase)
        α = makeArray(α_slice)
        MonRowInd = GetMonIndex(MonBase, α)
        for β_slice in eachrow(MonBase)
            β = makeArray(β_slice)
            MonColInd = GetMonIndex(MonBase, β)
            γ = α + β
            if haskey(MonMat, γ)
                # MonMat[tempEntry]  = push!(MonMat[γ], [α, β]) #  Index by exponents
                MonMat[γ]  = push!(MonMat[γ], (MonRowInd, MonColInd)) # Index by position in Monmial vector.

            else
                # MonMat[tempEntry]  = [[α, β]] #  Index by exponents
                MonMat[γ]  = [(MonRowInd, MonColInd)]
            end
        end
    end
    return MonMat
end
# @testset "GenMomMatDictComplex" begin
#     n = 2
#     t = 2
#     MonExp = GenMomMatDict(n,t*2)
#     MonBase = GenMon(n,t)
#     for key in keys(MonExp)
#         for val in MonExp[key]
#             @test   makeArray(key) == makeArray(MonBase[val[1],:] + MonBase[val[2],:])
#         end
#     end
# end
