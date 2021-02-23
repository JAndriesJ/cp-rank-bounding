using LinearAlgebra # For the matrix manipulations and preprocessing.
# using JuMP # For the optimization frame work.
include("moment_utils.jl")

# """This is just for testing puposes and sould not be used."""
# function make_dummy_var(n,t)
#     model        = Model(Mosek.Optimizer)
#     list_of_keys = make_mom_expo_keys(n, t)
#     @variable(model,Lx[list_of_keys] )
#     return Lx
# end


"""
input: A(Data matrix),t(integer),Lx(JuMP variable)
output: dictionary: keys:
comment: L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)
"""
function make_dag_con(A,t,Lx)
    n = size(A)[1]
    dag_con = Dict()

    # L(1⋅u) ≧ 0 for u ∈ [x]₂ₜ
    mom₂ₜ     = make_mon_expo(n,2*t)
    nb2t      = lastindex(mom₂ₜ)
    dag_con[(0,0)] = [ 1*Lx[mom₂ₜ[i] ] for i in  1:nb2t]


    deg_g = 2
    mom₂ₜ₋₂    = make_mon_expo(n,2*t-deg_g)
    nb2t₋₂ = lastindex(mom₂ₜ₋₂)
    @assert nb2t₋₂ > 1 "Contraints do not exist for t =1!"
    for k in 1:n
        eₖ = get_std_base_vec(n,k)
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        sqrMₖₖ = sqrt(A[k,k])
        dag_con[(k,k)] = [ sqrMₖₖ*Lx[eₖ + mom₂ₜ₋₂[i]] - Lx[2*eₖ + mom₂ₜ₋₂[i]] for i in  1:nb2t₋₂ ]
        for h in (k+1):n
            eₕ = get_std_base_vec(n,h)
            # Dagger constraints: L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
            dag_con[(k,h)] = [ A[k,h]*Lx[mom₂ₜ₋₂[i]] - Lx[eₖ + eₕ + mom₂ₜ₋₂[i]]     for i in  1:nb2t₋₂]
        end
    end
    return dag_con
end


"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (i,j) ∈ [n]²
                    vals:
comment: L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ) M_2
(M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function make_loc_con(A,t,Lx)
    n       = size(A)[1]
    LMB     = make_mon_expo(n, t - 1)
    nb_mon  = size(LMB)[1]
    loc_con = Dict()
    for k in 1:n
        eₖ = get_std_base_vec(n,k)
        # Constraint: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        sqrAₖₖ = sqrt(A[k,k])
        loc_con[(k,k)] = [sqrAₖₖ*Lx[LMB[i] + LMB[j] + eₖ] - Lx[LMB[i] + LMB[j] + 2*eₖ]   for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = get_std_base_vec(n,h)
            # Constraint: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            loc_con[(k,h)] = [ A[k,h]*Lx[LMB[i] + LMB[j]] - Lx[LMB[i] + LMB[j] + eₖ + eₕ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return loc_con
end


"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (h,k) ∈ [n]², h ≠ k
                    values:
"""
function make_xx_con(A,t,Lx)
    n       = size(A)[1]
    LMB     = make_mon_expo(n, t - 1)
    nb_mon  = size(LMB)[1]
    xx_con  = Dict()
    for k in 1:n
        eₖ = get_std_base_vec(n,k)
        for h in (k+1):n
            eₕ = get_std_base_vec(n,h)
            # Localizing xx constraint: M(xₖxₕ⋅L)
            xx_con[(k,h)] = [ Lx[LMB[i] + LMB[j] + eₖ + eₕ] for i in 1:nb_mon, j in 1:nb_mon ]
        end
    end
    return xx_con
end

## Tensor constraints

"""
Input: A(data matrix),t(Integer),Lx(JuMP variable)
Output: L((M-([x]₌₁[x]₌₁ᵀ))⊗([x]₌ₗ[x]₌ₗᵀ)))for l ∈ 0,1,...,t-1.
= M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))
"""
function make_weakG_con(A,t,Lx)
    n = size(A)[1]
    weakG_con = Dict()
    LMBexp_1 =  make_mon_expo_mat(n,1,false)
    for ℓ in 1:(t-1)
        LMBexp_ℓ          = make_mon_expo_mat(n,ℓ,false)   #exponents of [x]₌ₗ[x]₌ₗᵀ
        LMBexp_1ℓ         = var_kron(LMBexp_1,LMBexp_ℓ)    #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
        LMB_ℓ             = index_to_var(Lx,LMBexp_ℓ)      # L([x]₌ₗ[x]₌ₗᵀ)
        LMB_1ℓ            = index_to_var(Lx,LMBexp_1ℓ)     # L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))

        A_tens_LMB_ℓ      = kron(A,LMB_ℓ)                  # M⊗L([x]₌ₗ[x]₌ₗᵀ)
        weakG_con[ℓ]      = A_tens_LMB_ℓ - LMB_1ℓ          # M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ) ⪰ 0, ℓ ∈ 0,1,t-deg(g)/2
    end
    return weakG_con
end

"""M(G ⊗ L) ⪰ 0 constraints
input: A(data matrix),t(Integer),Lx(JuMP variable)
Assumption: G = A-[x]₌₁[x]₌₁ᵀ
output: A⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
"""
function make_G_con(A,t,Lx)
    n = size(A)[1]
    LMBexp_1           = make_mon_expo_mat(n,1,false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexpₜ₋₁          = make_mon_expo_mat(n,t-1,true)#exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁             = index_to_var(Lx,LMBexpₜ₋₁)    #L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp_1ₜ₋₁        = var_kron(LMBexp_1,LMBexpₜ₋₁)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB_1ₜ₋₁           = index_to_var(Lx,LMBexp_1ₜ₋₁)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB_1ₜ₋₁             # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
                                                  #L((A - ([x]₌₁[x]₌₁ᵀ))⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    #p = gen_tens_perm(n,t-1)
    return G_con#[p,p]
end

"""M(G ⊗ L) ⪰ 0 constraints"""
function make_G_con2(M,t,Lx)
    n = size(M)[1]

    LMBexp_00         = make_mon_expo_mat(n,0,false) #exponents of [x]₌₀[x]₌₀ᵀ
    LMBexp_10         = make_mon_expo_mat(n,(1,0),false) #exponents of [x]₌₁[x]₌₀ᵀ
    LMBexp_01         = make_mon_expo_mat(n,(0,1),false) #exponents of [x]₌₀[x]₌₁ᵀ
    LMBexp_11         = make_mon_expo_mat(n,1,false) #exponents of [x]₌₁[x]₌₁ᵀ

    B_exp_00          = var_kron(LMBexp_11,LMBexp_00 )   # ([x]₌₁[x]₌₁ᵀ)⊗([x]₌₀[x]₌₀ᵀ)
    B_exp_10          = var_kron(LMBexp_11,LMBexp_10)   # ([x]₌₁[x]₌₁ᵀ)⊗([x]₌₁[x]₌₀ᵀ)
    B_exp_01          = var_kron(LMBexp_11,LMBexp_01)   # ([x]₌₁[x]₌₁ᵀ)⊗([x]₌₀[x]₌₁ᵀ)
    B_exp_11          = var_kron(LMBexp_11,LMBexp_11)   # ([x]₌₁[x]₌₁ᵀ)⊗([x]₌₁[x]₌₁ᵀ)

    LMBexp            = vcat(hcat(B_exp_00,B_exp_01),hcat(B_exp_10,B_exp_11 ))
    B                 = index_to_var(Lx,LMBexp)

    A_00              = kron(M,index_to_var(Lx,LMBexp_00)) # M⊗L([x]₌₀[x]₌₀ᵀ)
    A_10              = kron(M,index_to_var(Lx,LMBexp_10)) # M⊗L([x]₌₁[x]₌₀ᵀ)
    A_01              = kron(M,index_to_var(Lx,LMBexp_01)) # M⊗L([x]₌₀[x]₌₁ᵀ)
    A_11              = kron(M,index_to_var(Lx,LMBexp_11)) # M⊗L([x]₌₁[x]₌₁ᵀ)

    A                 =vcat(hcat(A_00,A_01),hcat(A_10,A_11 ))
    GTensLCons = A - B
    return GTensLCons
end

"""
This function makes the permutation that links make_G_con and make_G_con2
    only holds for t = 2
"""
function gen_tens_perm(n)
    nar = [1 + j*(n+1) for j in 0:(n-1)]
    li = []
    for j in 1:(n*(n+1))
        if !(j in nar)
            append!(li,j)
        end
    end
    return [nar;li]
end

function bseq(n,k)
    li = Dict()
    li[0] = [1]
    for k_t in 0:k-1
        li[k_t+1] = [(1+binomial(n+k_t,k_t)):binomial(n+k_t+1,k_t+1)...]
    end
    return li
end


function gen_tens_perm(n,k)
    D = bseq(n,k)
    S = binomial(n+k,k)
    Li = []
    for key in keys(D)
        li = []
        for j in 0:(n-1)
            append!(li,D[key] .+ S*j)
        end
        Li = append!(Li,li)
    end
    return Li
end




## Utility
"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """
function index_to_var(var, index_array)
    sub = α -> var[α]
    var_array = sub.(index_array)
    return var_array
end



## Unused.

"""Not used anymore"""
function make_G_con1(LocConDict, M, LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    WeakGTensLConsDict = Dict()
    for k in 1:n
        eₖ = get_standard_base(n,k)
        # Diagonal blocks of the constraints:  M((Aₖₖ - xₖ²)⋅L)
        WeakGTensLConsDict[k] = [ M[k,k]*x[tra(LMB[i,:] + LMB[j,:])] - x[tra(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
    end


    GCon = []
    GConRow = []
    blockSize = size(LocConDict[1, 2])
    for i in 1:n
        for j in 1:n
            if j == 1
                if j == i
                    GConRow = WeakGTensLConsDict[1]
                else
                    GConRow = LocConDict[(min(i,1), max(i,1))]
                end
            else
                if  j == i
                    GConRow = hcat(GConRow, WeakGTensLConsDict[j])
                else
                    GConRow = hcat(GConRow, LocConDict[(min(i,j), max(i,j))])
                end
            end

        end
        # Just the first row
        if i == 1
            GCon = GConRow
        elseif i > 1
            GCon = vcat(GCon, GConRow)
        end
    end
    return GCon
end










"""Not used becuase it takes too long"""
function make_tens_constraint(M,t,Lx)
    LMBexp_1          = make_mon_expo_mat(n,1,false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexp_11         = var_kron(LMBexp_1,LMBexp_1)
    xxℓ               = var_self_kron(LMBexp_11,t)
    Lxxℓ              = index_to_var(Lx,xxℓ) # L( ([x]₌₁[x]₌₁ᵀ)^⊗ℓ) for ℓ = 0,1,..,t
    Mℓ                = self_kron(M,t)  # M^⊗ℓ for ℓ = 0,1,..,t)

    tens_constraint = Dict()
    for key in keys(Mℓ)
        tens_constraint[key] = Mℓ[key] - Lxxℓ[key]
    end
    return tens_constraint
end























# MakeGTensLConsMat(M,t,Lx) - MakeGTensLConsMat1(LocConDict, M, LMB,x)


"""This is never used """
function genCPCons(M,LMB,x)
    n = size(M)[1]
    LocConDict            = genCP_localizing_Constraints(M,LMB,x)
    DagConDict            = genCP_dagger_Constraints(M,LMB,x)
    LocConXXDict          = genCP_XX_Constraints(M,LMB,x)

    weakGTensLConsMat     = genCPweakGTensLCons(M,LMB,x)
    GTensLConsMat         = MakeGTensLConsMat(LocConDict, weakGTensLConsMat, n)

    return LocConDict, DagConDict, XXConDict, weakGTensLConsDict, GTensLConsMat
end

################################# MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, true)







################################# TO DO

"""
Tnesor constraints (L((ww')^c))w,w'∈<x>=l ≤  A⊗l for all integers 2 ≤ l ≤ t
"""
function GenTensConst(args)
    body
end













# Constraints matrices and vectors for cp rank:
# """
# The general form of a localizing matrix coming from CP rank is as follows:
# L ≥ 0 on M₂ₜ(S^cp_A )
# (M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ)
# M_2
# or
# (M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ)
#
# The dagger constraints are
# L(gu) ≥ 0 for g ∈ {1} ∪ S^cp_A and u ∈ [x]2t−deg(g)
#
# M(G ⊗ L) ⪰ 0 constraints
#
# """
# function genCPCons(M,LMB,x)
#     n = size(M)[1]
#     nb_mon  = size(LMB)[1]
#     LocConDict = Dict()
#     gConDict = Dict()
#     LocConXXDict = Dict()
#     DagConDict = Dict()
#     for k in 1:n
#         eₖ = LMB[(k+1),:]
#         # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
#         LocConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:] + LMB[j,:] + eₖ)] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
#         # Diagonal blocks of the constriants:  M((Aₖₖ - xₖ²)⋅L)
#         gConDict[k] = [M[k,k]*x[transpose(LMB[i,:] + LMB[j,:] )] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
#         # Localizing xx constriant: M(xₖ²⋅L)
#         LocConXXDict[(k,k)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ 2*eₖ ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
#         # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅m) ≧ 0 ??????????
#         DagConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:]  + LMB[(k+1),:] )] - x[transpose(LMB[i,:]  + 2*LMB[(k+1),:] )] for i in 1:nb_mon]
#         for h in (k+1):n
#             eₕ = LMB[(h+1),:]
#             # Constriant: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
#             LocConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
#             # Localizing xx constriant: M(xₖxₕ⋅L)
#             LocConXXDict[(k,h)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
#             # Dagger constraints: L((Aₖₕ  - xₖ²))  ??????????
#             DagConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
#         end
#     end
#     return LocConDict, DagConDict, LocConXXDict, gConDict
# end
