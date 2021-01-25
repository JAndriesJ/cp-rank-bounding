include("moment_utils.jl")

# Constraints matrices and vectors for cp rank:
"""
The general form of a localizing matrix coming from CP rank is as follows:
L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ)
M_2
or
(M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ)

The dagger constraints are
L(gu) ≥ 0 for g ∈ {1} ∪ S^cp_A and u ∈ [x]2t−deg(g)

M(G ⊗ L) ⪰ 0 constraints

"""
function genCPCons(M,LMB,x)
    n = size(M)[1]
    Loc_nb_mon  = size(LMB)[1]
    LocConDict = Dict()
    gConDict = Dict()
    LocConXXDict = Dict()
    DagConDict = Dict()
    for k in 1:n
        eₖ = LMB[(k+1),:]
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        LocConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:] + LMB[j,:] + eₖ)] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]
        # Diagonal blocks of the constriants:  M((Aₖₖ - xₖ²)⋅L)
        gConDict[k] = [M[k,k]*x[transpose(LMB[i,:] + LMB[j,:] )] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]
        # Localizing xx constriant: M(xₖ²⋅L)
        LocConXXDict[(k,k)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ 2*eₖ ) ] for i in 1:Loc_nb_mon,  j in 1:Loc_nb_mon ]
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅m) ≧ 0 ??????????
        DagConDict[(k,k)] = [sqrt(M[k,k])*x[transpose(LMB[i,:]  + LMB[(k+1),:] )] - x[transpose(LMB[i,:]  + 2*LMB[(k+1),:] )] for i in 1:Loc_nb_mon]
        for h in (k+1):n
            eₕ = LMB[(h+1),:]
            # Constriant: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            LocConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:Loc_nb_mon,  j in 1:Loc_nb_mon ]
            # Localizing xx constriant: M(xₖxₕ⋅L)
            LocConXXDict[(k,h)] = [ x[transpose(LMB[i,:] + LMB[j,:]+ eₖ + eₕ) ] for i in 1:Loc_nb_mon,  j in 1:Loc_nb_mon ]
            # Dagger constraints: L((Aₖₕ  - xₖ²))  ??????????
            DagConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:Loc_nb_mon,  j in 1:Loc_nb_mon ]
        end
    end
    return LocConDict, DagConDict, gConDict , LocConXXDict
end

""" M(G ⊗ L) ⪰ 0 """
function MakeG(LocConDict, GConDict,n)
    GCon = []
    GConRow = []
    for i in 1:n
        for j in 1:n
            if j == 1
                if j == i
                    GConRow = GConDict[1]
                else
                    GConRow = LocConDict[(min(i,1), max(i,1))]
                end
            else
                if  j == i
                    GConRow = hcat(GConRow, GConDict[j])
                else
                    GConRow = hcat(GConRow, LocConDict[(min(i,j), max(i,j))])
                end
            end
        end
        # println(size(GConRow))
        if i == 1
            GCon = GConRow
        elseif i > 1
            GCon = vcat(GCon, GConRow)
        end
    end
    return GCon
end
#@testset "genCPCons" begin


"""
Tnesor constraints (L((ww')^c))w,w'∈<x>=l ≤  A⊗l for all integers 2 ≤ l ≤ t
"""
function GenTensConst(args)
    body
end
