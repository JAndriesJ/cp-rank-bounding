include("moment_utils.jl")

""" sep rank constraints """
function genSEPCons(ρ,LMB,x)
    d = sqrt(size(ρ)[1])
    √ρₘₐₓ = sqrt(max(diag(ρ)))
    √trρ = sqrt(tr(ρ))
    Loc_nb_mon  = size(LMB)[1]
    LocConXDict = Dict()
    LocConYDict = Dict()
    DagConDict = Dict()
    gConDict = Dict()
    for k in 1:n
        eₖ = LMB[(k+1),:]
        eₖ₊ₙ = LMB[(k+1+n),:]
        eₖ₊₂ₙ  = LMB[(k+1+2*n),:]
        eₖ₊₃ₙ = LMB[(k+1+3*n),:]
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)
        # Sᶜᵖ := {√ρₘₐₓ -  xᵢ̄xⱼ} ∪ {√ρₘₐₓ -  yᵢ̄yⱼ}
            # polynomials in {√ρₘₐₓ -  xᵢ̄xⱼ}
        LocConXDict[(k,k)] = [√ρₘₐₓ- x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₖ₊ₙ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]
            # polynomials in {√ρₘₐₓ -  yᵢ̄yⱼ}
        LocConYDict[(k,k)] = [√ρₘₐₓ- x[transpose(LMB[i,:] + LMB[j,:] + eₖ₊₂ₙ + eₖ₊₃ₙ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]


        # Sᶜᵖ := {√trρ -  xx̄} ∪ {√trρ -  yȳ}
        # polynomials in {√ρₘₐₓ -  xᵢ̄xⱼ  }
            #LocConXDict[(k,k)] = [√trρ- x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₖ₊ₙ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]
        # polynomials in {√ρₘₐₓ -  yᵢ̄yⱼ  }
            #LocConYDict[(k,k)] = [√trρ- x[transpose(LMB[i,:] + LMB[j,:] + eₖ₊₂ₙ + eₖ₊₃ₙ)] for i in 1:Loc_nb_mon, j in 1:Loc_nb_mon ]

        for h in (k+1):n
            eₕ = LMB[(h+1),:]

            # Localizing g constriant : M((ρᵢᵢ⋅ⱼⱼ⋅ - xᵢ̄xⱼyᵢ̄yⱼ⋅)⋅L)
            LocConDict[(k,h)] = [ M[k,h]*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ + eₕ) ] for i in 1:Loc_nb_mon,  j in 1:Loc_nb_mon ]

        end
    end
    return LocConDict, DagConDict, gConDict , LocConXXDict
end


""" M(G ⊗ L) ⪰ 0 """
function MakeGρ(LocConDict, GConDict,n)
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







"""
Tnesor constraints (L((ww')^c))w,w'∈<x>=l ≤  A⊗l for all integers 2 ≤ l ≤ t
"""
function GenTensConst(args)
    body
end


# using LinearAlgebra # For the matrix manipulations and preprocessing.
# using JuMP # For the optimization frame work.
# using MosekTools # The solver that we use.
