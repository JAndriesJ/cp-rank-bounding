using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.

include("moment_utils.jl")

"""L(gu) ≥ 0 for g ∈ {1} ∪ S^cp_A and u ∈ [x]2t−deg(g)"""
function genCP_dagger_Constraints(M,t,Lx)
    n = size(M)[1]
    DagConDict = Dict()

    # L(1⋅u) ≧ 0 for u ∈ [x]₂ₜ
    LMB_full = make_mon_expo(n,2*t)
    nb_mon_full  = size(LMB_full)[1]
    DagConDict[(0,0)] = [ 1*Lx[tra(LMB_full[i,:]) ]  for i in  1:nb_mon_full]

    deg_g = 2
    LMB = make_mon_expo(n, 2*t - deg_g)
    nb_mon = size(LMB)[1]
    @assert nb_mon > 1 "Contraints do not exist for t =1!"
    for k in 1:n
        eₖ = get_standard_base(n,k)
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        DagConDict[(k,k)] = [ sqrt(M[k,k])*Lx[tra(eₖ + LMB[i,:]) ] - Lx[tra(2*eₖ + LMB[i,:])] for i in  1:nb_mon ]

        for h in (k+1):n
            eₕ = get_standard_base(n,h)
            # Dagger constraints: L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
            DagConDict[(k,h)] = [ M[k,h]*Lx[tra(LMB[i,:] )] - Lx[tra(eₖ + eₕ + LMB[i,:])]     for i in  1:nb_mon]
            # MISTAKE below
            # DagConDict[(k,h)] = [ M[k,h]*x[tr(LMB[i,:] + LMB[j,:]  )] - x[tr(eₖ + eₕ + LMB[i,:] + LMB[j,:] )]     for i in  1:nb_mon, j in  1:nb_mon  ]
        end
    end
    return DagConDict
end

""" L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ) M_2 or (M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function genCP_localizing_Constraints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    LocConDict = Dict()
    for k in 1:n
        eₖ = get_standard_base(n,k)
        # Constraint: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        LocConDict[(k,k)] = [sqrt(M[k,k])*x[tra(LMB[i,:] + LMB[j,:] + eₖ)] - x[tra(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = get_standard_base(n,h)
            # Constraint: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            LocConDict[(k,h)] = [ M[k,h]*x[tra(LMB[i,:] + LMB[j,:])] - x[tra(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return LocConDict
end

function genCP_XX_Constraints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    XXConDict = Dict()
    for k in 1:n
        eₖ = get_standard_base(n,k)
        # Localizing xx constraint: M(xₖ²⋅L)
        XXConDict[(k,k)] = [ x[tra(LMB[i,:] + LMB[j,:] + 2*eₖ ) ] for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = get_standard_base(n,h)
            # Localizing xx constraint: M(xₖxₕ⋅L)
            XXConDict[(k,h)] = [ x[tra(LMB[i,:] + LMB[j,:] + eₖ + eₕ ) ] for i in 1:nb_mon, j in 1:nb_mon ]
        end
    end
    return XXConDict
end

## Tensor constraints

"""Takes an array of exponents α's and gives array of same shape L(xᵅ)  """
function index_to_var(var, index_array)
    n,m = size(index_array)
    var_array = Array{Any}(undef, n,m)
    for i in 1:n
        for j in 1:m
            var_array[i,j] = var[tra(index_array[i,j])]
        end
    end
    return var_array
end

# function assemble_dict(dict_of_blocks)
#     n = Int(sqrt(length(keys(dict_of_blocks))))
#     # check if this is integer.
#     if n == 1
#         return dict_of_blocks[1,1]
#     end
#     block = []
#     row_block = []
#     for i in 1:n
#         for j in 1:n
#             if j == 1
#                 row_block = dict_of_blocks[i, j]
#                 # if j == i
#                 #     row_block = dict_of_blocks[1,1]
#                 # else
#
#                 # end
#             else
#                 row_block = hcat(row_block, dict_of_blocks[i,j])
#             end
#
#         end
#         # Just the first row
#         if i == 1
#             block = row_block
#         elseif i > 1
#             block = vcat(block, row_block)
#         end
#     end
#     return block
# end
"""
Input: M,t,x
Output: L((M-([x]₌₁[x]₌₁ᵀ))⊗([x]₌ₗ[x]₌ₗᵀ)))
= M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))
for l ∈ 0,1,...,t-1.
"""
function genCPweakGTensLCons(M,t,Lx)
    n = size(M)[1]
    weakGTensLCons = Dict()
    LMBexp_1 =  make_mon_expo_mat(n,1,false)
    for ℓ in 0:(t-1)
        LMBexp_ℓ          = make_mon_expo_mat(n,ℓ,false)   #exponents of [x]₌ₗ[x]₌ₗᵀ
        LMBexp_1ℓ         = var_kron(LMBexp_1,LMBexp_ℓ)    #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
        LMB_ℓ             = index_to_var(Lx,LMBexp_ℓ)      # L([x]₌ₗ[x]₌ₗᵀ)
        LMB_1ℓ            = index_to_var(Lx,LMBexp_1ℓ)     # L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))

        M_tens_LMB_ℓ      = kron(M,LMB_ℓ)                  # M⊗L([x]₌ₗ[x]₌ₗᵀ)
        weakGTensLCons[ℓ] = M_tens_LMB_ℓ - LMB_1ℓ         # M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
    end
    return weakGTensLCons
end
"""M(G ⊗ L) ⪰ 0 constraints"""
function MakeGTensLConsMat(M,t,Lx)
    n = size(M)[1]
    LMBexp_1          = make_mon_expo_mat(n,1,false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexp_t          = make_mon_expo_mat(n,1,true)  #exponents of [x]ₜ[x]ₜᵀ
    LMB_t             = index_to_var(Lx,LMBexp_t)    #L([x]ₜ[x]ₜᵀ)

    LMBexp_1t         = var_kron(LMBexp_1,LMBexp_t)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
    LMB_1t            = index_to_var(Lx,LMBexp_1t)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ[x]ₜᵀ))

    GTensLCons = kron(M,LMB_t) - LMB_1t              # M⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ[x]ₜᵀ)
    return GTensLCons
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
"""Not used"""
function MakeGTensLConsMat1(LocConDict, M, LMB,x)
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
