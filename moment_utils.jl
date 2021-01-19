#Do you have a moment? (Goal: implement xi_cp of the FOCM paper. Eventually add the G-constraints)
using LinearAlgebra # For the matrix manipulations and preprocessing.
using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
import .mat_repo
include("..\\matrix-menagerie\\./mat_repo.jl")

# abbreviations for convenience
tra = transpose
cl = clearconsole

"""
input: n (integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of integers)
coments: Has a quadratic loss in runtime because for degree t we compute degree 1 t times...
"""
function make_mon_expo_arr(n::Int,t::Int, isLeq::Bool = true)
    if t == 0 # [x]₌₀
        return zeros(Int32,1,n)
    else # [x]₌ₜ
        temp = make_mon_expo_arr(n,t-1,isLeq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if isLeq # [x]≦ₜ
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end
# println(make_mon_expo_arr(3,2,false))

"""
input: n (integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of arrays of integers)
"""
function make_mon_expo(n::Int,t::Int, isLeq::Bool = true)
    mon_expo_arr = make_mon_expo_arr(n,t,isLeq)
    mon_expo     = [r  for r in  eachrow(mon_expo_arr)]
    return mon_expo
end
# println(make_mon_expo(3,2,false))


"""
input: n (integer),k(integer)
output: eₖ ∈ {0,1}ⁿ i.e. the standard basis vector
"""
function get_std_base_vec(n::Int,k::Int)
    mon_expo = make_mon_expo(n,1,false)
    eₖ = mon_expo[k]
    return eₖ
end
# println(get_std_base_vec(3,3))


"""
input: n(integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ[x]≦ₜᵀ or [x]₌ₜ[x]₌ₜᵀ where , x = (x₁,x₂,...,xₙ)
"""
function make_mon_expo_mat(n::Int,t::Tuple{Int64,Int64},isLeq::Bool = true)
    mon1      = make_mon_expo(n,t[1], isLeq)
    mon2      = make_mon_expo(n,t[2], isLeq)
    nb_mon1   = length(mon1)
    nb_mon2   = length(mon2)
    xxᵀₜ_vec = [ mon1[i]  + mon2[j] for i in 1:nb_mon1 for j in 1:nb_mon2]

    xxᵀₜ     = reshape(xxᵀₜ_vec, (nb_mon1, nb_mon2) )
    return xxᵀₜ
end

function make_mon_expo_mat(n::Int,t::Int,isLeq::Bool = true)
    xxᵀₜ     = make_mon_expo_mat(n,(t,t),isLeq)
    return xxᵀₜ
end

# println(make_mon_expo_mat(3,(2,1),false))



"""
input: B ∈ (N_t^n)ⁿᵐ (multi-index array), α ∈ N_t^n (multi-index)
output:  row index of B containing α or error message
comment: For some reason the i-th row of an array is returned as a column vector (the transpose remedies this)
"""
function get_mon_index(B,α)
    try
        return findall(Bool[tra(B[i]) == α for i =1:size(B,1)])[1]
    catch
        println("There is no such entry.")
    end
end
# mom = make_mon_expo(3,2,false)
# println(get_mon_index(mom,[1 1 0]))


"""
input:  n(integer),t(integer)
output: dictionary: keys: unique exponents in [x]≦ₜ[x]≦ₜᵀ
                    values: indeces in [x]≦ₜ[x]≦ₜᵀ corresponding to key as exponent
comment: Dictionary: keys γ ∈ N_2t^n, values are indeces in Moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ
"""
function make_mom_expo_mat_dict(n::Int,t::Int)
    mon_vec = make_mon_expo(n,t)
    mom_mat_dict = Dict()

    for α in mon_vec
        row_index = get_mon_index(mon_vec,tra(α))
        for β in mon_vec
            col_index = get_mon_index(mon_vec, tra(β))
            γ = α + β
            if haskey(mom_mat_dict, γ)
                mom_mat_dict[γ]  = push!(mom_mat_dict[γ], (row_index, col_index)) # Index by position in Monmial vector.
            else
                mom_mat_dict[γ]  = [(row_index, col_index)]
            end
        end
    end
    return mom_mat_dict
end

## Utility

"""
input: Dictionary:
    keys: i,j ∈ [n], n::Int
    values: square array A_{i,j}
output: Array A = (A_{i,j})_i,j
"""
function assemble_dict(dict_of_blocks)
    n = Int(sqrt(length(keys(dict_of_blocks))))
    if n == 1
        return dict_of_blocks[1,1]
    end
    block = []
    row_block = []
    for i in 1:n
        for j in 1:n
            if j == 1
                row_block = dict_of_blocks[i, j]
            else
                row_block = hcat(row_block, dict_of_blocks[i,j])
            end
        end

        if i == 1
            block = row_block
        elseif i > 1
            block = vcat(block, row_block)
        end

    end
    return block
end

##
"""
input: A,B (arrays of integer tupples)
output: exponent array of A ⊗ B
"""
function var_kron(A,B)
    n1,n2 = size(A)

    D = Dict()
    for i in 1:n1
        for j in 1:n2
            C = repeat( [A[i,j]] , inner = (1,1), outer = size(B))
            D[i,j] = C + B
        end
    end
    return assemble_dict(D)
end

"""
input: A(array of integer tupples), ℓ(Integer)
output: exponent array of A ⊗ ...⊗ A (ℓ-times)
"""
function var_self_kron(A,ℓ)
    A_tens = Dict()
    A_tens[0] = A
    for i in 1:ℓ
        A_tens[i] = var_kron(A,A_tens[i-1])
    end
    return A_tens
end
