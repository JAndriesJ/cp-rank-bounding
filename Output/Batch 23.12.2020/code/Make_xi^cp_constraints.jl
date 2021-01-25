using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.

include("moment_utils.jl")

""" L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ) M_2 or (M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function genCP_localizing_Constriaints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    LocConDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        LocConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:] + LMB[j,:] + eₖ)] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Constriant: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            LocConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return LocConDict
end

"""L(gu) ≥ 0 for g ∈ {1} ∪ S^cp_A and u ∈ [x]2t−deg(g)"""
function genCP_dagger_Constriaints(M,LMB,x)
    @assert size(LMB)[1] > 1 "Contraints do not exist for t =1!"
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    DagConDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅m) ≧ 0
        DagConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:]  + LMB[(k+1),:] )] - x[transpose(LMB[i,:]  + 2*LMB[(k+1),:] )] for i in 1:nb_mon]
        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Dagger constraints: L((Aₖₕ  - xₖ²)) ≧ 0 
            DagConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return DagConDict
end

function genCP_XX_Constriaints(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    XXConDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Localizing xx constriant: M(xₖ²⋅L)
        XXConDict[(k,k)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ 2*eₖ ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = standardBase(n,h)
            # Localizing xx constriant: M(xₖxₕ⋅L)
            XXConDict[(k,h)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ eₖ + eₕ) ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return XXConDict
end

function genCPpreWeakGTensLCons(M,LMB,x)
    n = size(M)[1]
    nb_mon  = size(LMB)[1]
    weakGTensLConsDict = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Diagonal blocks of the constriants:  M((Aₖₖ - xₖ²)⋅L)
        weakGTensLConsDict[k] = [M[k,k]*x[transpose(LMB[i,:] + LMB[j,:] )] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
    end
    return weakGTensLConsDict
end

"""M(G ⊗ L) ⪰ 0 constraints"""
function MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, isWeak=false)

    GCon = []
    GConRow = []
    blockSize = size(LocConDict[1, 2])
    for i in 1:n
        for j in 1:n
            if j == 1
                if j == i
                    GConRow = preWeakGTensLConsDict[1]
                else
                    if isWeak
                        GConRow = zeros(blockSize)
                    else
                        GConRow = LocConDict[(min(i,1), max(i,1))]
                    end
                end
            else
                if  j == i
                    GConRow = hcat(GConRow, preWeakGTensLConsDict[j])
                else
                    if isWeak
                        GConRow = hcat(GConRow, zeros(blockSize))
                    else
                        GConRow = hcat(GConRow, LocConDict[(min(i,j), max(i,j))])
                    end
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


function genCPCons(M,LMB,x)
    n = size(M)[1]
    LocConDict            = genCP_localizing_Constriaints(M,LMB,x)
    DagConDict            = genCP_dagger_Constriaints(M,LMB,x)
    LocConXXDict          = genCP_XX_Constriaints(M,LMB,x)
    preWeakGTensLConsDict = genCPpreWeakGTensLCons(M,LMB,x)

    weakGTensLConsMat     = MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, true)
    GTensLConsMat         = MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, false)

    return LocConDict, DagConDict, XXConDict, weakGTensLConsDict, GTensLConsMat
end

#################################







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
