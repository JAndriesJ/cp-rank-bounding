#Do you have a moment? (Goal: implement xi_cp of the FOCM paper. Eventually add the G-constraints)
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

""" The name is self explanitory """
function standardBase(n,k)
    LMB = GenMon(n, 1)
    eₖ = LMB[(k+1),:]
    return eₖ
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



"""Convert a sub-array pointer to an array."""
function makeArray(thing)
    α = copy(thing)
    n = length(α)
    α = reshape(α,1,n)
    return α
end


"""Dictionary: keys γ ∈ N_2t^n, values are indeces in Moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ"""
function GenMomMatDict(n,t)
    MonBase = GenMon(n,t)
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


# """Dictionary: keys γ ∈ N_2t^n, values are indeces in Moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ"""
# function GenMomMatDictComplex(n,t)
#     @assert iseven(dt) "$t cannot be an odd number!"
#     MonBase = GenMon(n,t)
#     MonMat = Dict()
#     for α_slice in eachrow(MonBase)
#         α = makeArray(α_slice)
#         MonRowInd = GetMonIndex(MonBase, α)
#         for β_slice in eachrow(MonBase)
#             β = makeArray(β_slice)
#             MonColInd = GetMonIndex(MonBase, β)
#             γ = α + β
#             if haskey(MonMat, γ)
#                 # MonMat[tempEntry]  = push!(MonMat[γ], [α, β]) #  Index by exponents
#                 MonMat[γ]  = push!(MonMat[γ], (MonRowInd, MonColInd)) # Index by position in Monmial vector.
#
#             else
#                 # MonMat[tempEntry]  = [[α, β]] #  Index by exponents
#                 MonMat[γ]  = [(MonRowInd, MonColInd)]
#             end
#         end
#     end
#     return MonMat
# end
