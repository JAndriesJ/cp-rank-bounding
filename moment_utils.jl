#Do you have a moment? (Goal: implement xi_cp of the FOCM paper. Eventually add the G-constraints)
using LinearAlgebra # For the matrix manipulations and preprocessing.
tra = transpose
cl = clearconsole

"""Input: n,t output: exponents of [x]≦ₜ."""
function make_mon_expo(n,t, isLeq = true)
    if t == 0
        return zeros(Int32,1,n)
    else
        temp = make_mon_expo(n,t-1,isLeq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if isLeq # For monomials of degree less than t
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end #Has a quadratic loss in runtime because for degree t we compute degree 1 t times...


"""Input: n,t output: exponents of [x]≦ₜ[x]≦ₜᵀ where , x = (x₁,x₂,...,xₙ)."""
function make_mon_expo_mat(n,t,isLeq = true)
    if typeof(t) == Int64
        t = (t,t)
    end
    #     mon      = make_mon_expo(n,t, isLeq)
    #     nb_mon   = size(mon)[1]
    #     xxᵀₜ_vec = [ mon[i,:]  + mon[j,:] for i in 1:nb_mon for j in 1:nb_mon]
    #     xxᵀₜ     = reshape(xxᵀₜ_vec, (nb_mon, nb_mon) )
    # elseif  t == Tuple{Int64,Int64}
    mon1      = make_mon_expo(n,t[1], isLeq)
    mon2      = make_mon_expo(n,t[2], isLeq)
    nb_mon1   = size(mon1)[1]
    nb_mon2   = size(mon2)[1]
    xxᵀₜ_vec = [ mon1[i,:]  + mon2[j,:] for i in 1:nb_mon1 for j in 1:nb_mon2]
    xxᵀₜ     = reshape(xxᵀₜ_vec, (nb_mon1, nb_mon2) )
    # else
    #     @warn "In correct level specification."
    # end
    return xxᵀₜ
end
# @test make_mon_expo_mat(7,(1,1),true) == make_mon_expo_mat(7,1,true)


"""Input is a multi-index array N_t^n and α ∈ N_t^n, output is the row index of the
 monomial with power equal to the multi-index."""
function get_mon_index(M,α)
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
function make_mom__expo_mat_dict(n,t)
    MonBase = make_mon_expo(n,t)
    MonMat = Dict()
    for α_slice in eachrow(MonBase)
        α = makeArray(α_slice)
        MonRowInd = get_mon_index(MonBase, α)
        for β_slice in eachrow(MonBase)
            β = makeArray(β_slice)
            MonColInd = get_mon_index(MonBase, β)
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


""" The name is self explanitory """
function get_standard_base(n,k)
    LMB = make_mon_expo(n, 1)
    eₖ = LMB[(k+1),:]
    return eₖ
end


"""Input: Dictionary:
keys: i,j
values: square arrays A_{i,j}
Output: Array A = (A_{i,j})_i,j
"""
function assemble_dict(dict_of_blocks)
    n = Int(sqrt(length(keys(dict_of_blocks))))
    # check if this is integer.
    if n == 1
        return dict_of_blocks[1,1]
    end
    block = []
    row_block = []
    for i in 1:n
        for j in 1:n
            if j == 1
                row_block = dict_of_blocks[i, j]
                # if j == i
                #     row_block = dict_of_blocks[1,1]
                # else

                # end
            else
                row_block = hcat(row_block, dict_of_blocks[i,j])
            end

        end
        # Just the first row
        if i == 1
            block = row_block
        elseif i > 1
            block = vcat(block, row_block)
        end
    end
    return block
end

##
function var_kron(A,B)
    n1,n2 = size(A)
    m1,m2 = size(B)

    D = Dict()
    for i in 1:n1
        for j in 1:n2
            C = repeat( [A[i,j]] , inner = (1,1), outer = (m1,m2))
            D[i,j] = B + C
        end
    end
    return assemble_dict(D)
end

function var_self_kron(A,ℓ)
    A_tens = Dict()
    A_tens[0] = A
    for i in 1:ℓ
        A_tens[i] = var_kron(A,A_tens[i-1])
    end
    return A_tens
end


# A = [ [[1 0 0 2 0 1 0]]  [[1 2 0 0 0 0 0]]  [[0 1 0 2 0 0 1]]  [[0 0 0 1 2 0 0]] ]
# A = reshape(A,2,2)
# B = [ [[1 0 0 2 0 1 0]]  [[1 0 0 2 0 1 0]]  [[1 0 0 2 0 1 0]]  [[1 0 0 2 0 1 0]] ]
# B = reshape(B,2,2)
# D = var_kron(A,B)


function gen_mon_leq_mat(n,t)
    B = gen_mon_mat(n,t-2,true)

    g_mon_eq_mat = Dict()
    for i in 1:n
        eᵢ = standardBase(n,i)
        for j in 1:n
            eⱼ = standardBase(n,j)
            B = repeat( tra(eᵢ + eⱼ) , inner = (1,1), outer = (K,K))

            g_mon_eq_mat[(i,j)] = M[i,j]*mon_eq_mat - index_to_var(var, mon_eq_mat_ex + B)
        end
    end
end





""" """
function gen_g_mon_eq_mat(M,mon_eq_mat,mon_eq_mat_ex,var)
    n = size(M)[1]
    K = size(mon_eq_mat)[1]

    g_mon_eq_mat = Dict()
    for i in 1:n
        eᵢ = standardBase(n,i)
        for j in 1:n
            eⱼ = standardBase(n,j)
            B = repeat( tra(eᵢ + eⱼ) , inner = (1,1), outer = (K,K))

            g_mon_eq_mat[(i,j)] = M[i,j]*mon_eq_mat - index_to_var(var, mon_eq_mat_ex + B)
        end
    end
    return g_mon_eq_mat
end
