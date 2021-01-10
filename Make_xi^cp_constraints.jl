using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.

include("moment_utils.jl")



"""L(gu) ≥ 0 for g ∈ {1} ∪ S^cp_A and u ∈ [x]2t−deg(g)"""
function genCP_dagger_Constriaints(M,t,x)
    n = size(M)[1]

    deg_g = 2
    LMB = GenMon(n, 2*t - deg_g)
    @assert size(LMB)[1] > 1 "Contraints do not exist for t =1!"

    nb_mon  = size(LMB)[1]
    DagConDict = Dict()

    # L(1⋅u) ≧ 0 for u ∈ [x]₂ₜ
    LMB_full = GenMon(n, 2*t )
    nb_mon_full  = size(LMB_full)[1]
    DagConDict[(0,0)] = [ 1*x[tra(LMB_full[i,:]) ]  for i in  1:nb_mon_full]
    for k in 1:n
        eₖ = standardBase(n,k)
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        DagConDict[(k,k)] = [ sqrt(M[k,k])*x[tra(eₖ + LMB[i,:]) ] - x[tra(2*eₖ + LMB[i,:])] for i in  1:nb_mon ]

        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Dagger constraints: L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
            #g = ( M[k,h] - x[tr(eₖ + eₕ)]  )

            DagConDict[(k,h)] = [ M[k,h]*x[tra(LMB[i,:] )] - x[tra(eₖ + eₕ + LMB[i,:])]     for i in  1:nb_mon]
            # MISTAKE below
            # DagConDict[(k,h)] = [ M[k,h]*x[tr(LMB[i,:] + LMB[j,:]  )] - x[tr(eₖ + eₕ + LMB[i,:] + LMB[j,:] )]     for i in  1:nb_mon, j in  1:nb_mon  ]
        end
    end
    return DagConDict
end


""" L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ) M_2 or (M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function genCP_localizing_Constriaints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    LocConDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        LocConDict[(k,k)] = [sqrt(M[k,k])*x[tra(LMB[i,:] + LMB[j,:] + eₖ)] - x[tra(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Constriant: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            LocConDict[(k,h)] = [ M[k,h]*x[tra(LMB[i,:] + LMB[j,:])] - x[tra(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return LocConDict
end



function genCP_XX_Constriaints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    XXConDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Localizing xx constriant: M(xₖ²⋅L)
        XXConDict[(k,k)] = [ x[tra(LMB[i,:] + LMB[j,:] + 2*eₖ ) ] for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Localizing xx constriant: M(xₖxₕ⋅L)
            XXConDict[(k,h)] = [ x[tra(LMB[i,:] + LMB[j,:] + eₖ + eₕ ) ] for i in 1:nb_mon, j in 1:nb_mon ]
        end
    end
    return XXConDict
end


function genCPweakGTensLCons(M,t,x)
    n = size(M)[1]
    weakGTensLCons = Dict()
    for t_temp = 0:t/2
        mon_eq_mat_ex  = gen_mon_eq_mat(n,t_temp)
        mon_eq_mat     = index_to_var(x, mon_eq_mat_ex)
        g_mon_eq_mat   = gen_g_mon_eq_mat(M,mon_eq_mat,mon_eq_mat_ex,x)
        weakGTensLCons[t_temp] = g_mon_eq_mat
    end
    return weakGTensLCons
end




"""M(G ⊗ L) ⪰ 0 constraints"""
function MakeGTensLConsMat(LocConDict, M, LMB)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    WeakGTensLConsDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Diagonal blocks of the constriants:  M((Aₖₖ - xₖ²)⋅L)
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
    LocConDict            = genCP_localizing_Constriaints(M,LMB,x)
    DagConDict            = genCP_dagger_Constriaints(M,LMB,x)
    LocConXXDict          = genCP_XX_Constriaints(M,LMB,x)

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
