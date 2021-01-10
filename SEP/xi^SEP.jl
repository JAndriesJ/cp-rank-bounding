
using LinearAlgebra # For the matrix manipulations and preprocessing.
using JuMP # For the optimization frame work.
using Test # Utility module for testing if code works as intended.
using MosekTools # The solver that we use.

include("Matrix_Examples_SEP.jl") # Import matrices Examples we are interested in.
include("moment_utils.jl")


testCode = true
if testCode
    d₁ = 3
    d₂ = 3
    t = 3
    r = 3
    n =d₁ + d₂
    ρ = genRandGSEPmatrix(d₁, d₂, r)
    MomMatExp  = GenMomMatDict(n, t)
    MonBase = GenMon(n, t)
    model = Model(Mosek.Optimizer)
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )
end




""" L ≥ 0 on M₂ₜ(S^SEP_ρ )
M((√ρₘₐₓ - xₖ²)⋅L) ≥ 0
or
M((√ρₘₐₓ - yₖ²)⋅L) ≥ 0 """
function genSEP_localizing_Constriaints(ρ,t,x)
    sqrtρₘₐₓ = sqrt(maximum(diag(ρ)))
    n = Int(2*sqrt(size(ρ)[1]))
    LMB = GenMon(n, t-1)
    nb_mon  = size(LMB)[1]
    LocConDictX = Dict()
    LocConDictY = Dict()
    for k in 1:n
        eₖ = standardBase(n,k)
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√ρₘₐₓ - xₖ²)⋅L)
        LocConDictX[k] = [sqrtρₘₐₓ*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
        # Constriant: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√ρₘₐₓ - yₖ²)⋅L)
        LocConDictY[k] = [sqrtρₘₐₓ*x[transpose(LMB[i,:] + LMB[j,:])] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ)] for i in 1:nb_mon, j in 1:nb_mon ]
    end
    return LocConDictX, LocConDictY
end
if testCode
    @testset "genSEP_localizing_Constriaints" begin
        LocConDictX, LocConDictY = genSEP_localizing_Constriaints(ρ,t,x)
    end
end




""" L ≥ 0 on M₂ₜ(S^SEP_ρ )
 M((√Trρ - ||x||²)⋅L) ≥ 0
or
M((√Trρ - ||y||²)⋅L) ≥ 0 """
function genSEP_Tr_localizing_Constriaints(ρ,t,x)
    sqrtTrρ = sqrt(sum(diag(ρ)))
    n = Int(2*sqrt(size(ρ)[1]))
    LMB = GenMon(n, t-1)
    d₁,d₂ = Int(n/2),Int(n/2)
    nb_mon  = size(LMB)[1]

    LocConDictX = Dict()
    LocConDictY = Dict()

    normXsqr =  hcat(2*ones(1,d₁),zeros(1,d₂))
    normYsqr =  hcat(zeros(1,d₁),2*ones(1,d₂))

    normXsqr =  round.(Int, normXsqr)
    normYsqr =  round.(Int, normYsqr)

    # Constriant: L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Trρ - ||x||²)⋅L)
    LocConDictX = [sqrtTrρ*x[transpose(LMB[i,:] + LMB[j,:])] -  sum( [x[transpose(LMB[i,:] + LMB[j,:] + 2*standardBase(n,k) )] for k in 1:d₁ ])   for i in 1:nb_mon, j in 1:nb_mon]
    # Constriant:  L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Trρ - ||y||²)⋅L)
    LocConDictY = [sqrtTrρ*x[transpose(LMB[i,:] + LMB[j,:])] -  sum( [x[transpose(LMB[i,:] + LMB[j,:] + 2*standardBase(n,d₁ + k) )] for k in 1:d₂ ])    for i in 1:nb_mon, j in 1:nb_mon]

    return LocConDictX, LocConDictY
end
if testCode
    @testset "genSEP_Tr_localizing_Constriaints" begin
        LocConDictX, LocConDictY = genSEP_Tr_localizing_Constriaints(ρ,t,x)
    end
end

function genSEPGTensLConsDict(ρ,t,x,isWeak=false)
    @assert t > 1 "Hierarchy level too low. These constriants only hold for hiegher level steps of the hierarchy!"
    n = Int(2*sqrt(size(ρ)[1]))
    LMB = GenMon(n, t-2)
    d₁,d₂ = Int(n/2),Int(n/2)
    nb_mon  = size(LMB)[1]
    GTensLConsDict = Dict()
    # k₁,k₂,h₁,h₂ = rand(1:d₁,(1,4))
    # eₖ₁ = standardBase(n,k₁)
    # eₖ₂ = standardBase(n, k₂)
    # eₕ₁ = standardBase(n,d₁ + h₁)
    # eₕ₂ = standardBase(n,d₁ + h₂)
    for k₁ in 1:d₁
        eₖ₁ = standardBase(n,k₁)
        for k₂ in 1:d₁
            eₖ₂ = standardBase(n, k₂)
            for h₁ in 1:d₂
                eₕ₁ = standardBase(n,d₁ + h₁)
                for h₂ in 1:d₂
                    eₕ₂ = standardBase(n,d₁ + h₂)
                    if  isWeak
                        # Constriants:  M((ρₖₖₕₕ - xₖ²yₕ²)⋅L) on Diagonal blocks...
                        if k₁ == k₂ && h₁ == h₂
                            GTensLConsDict[k₁,k₂,h₁,h₂] =  [ρ[d₁*(k₁-1) + k₂, d₂*(h₁-1) + h₂]*x[transpose(LMB[i,:] + LMB[j,:] )] - x[transpose(LMB[i,:] + LMB[j,:] + 2*eₖ₁ + 2*eₕ₁)] for i in 1:nb_mon, j in 1:nb_mon ]
                        else
                            #... and Zero otherwise
                            GTensLConsDict[k₁,k₂,h₁,h₂] = zeros(nb_mon,nb_mon)
                        end
                    else
                        # Constriants:  M((ρₖ₁ₖ₂ₕ₁ₕ₂ - xₖ₁xₖ₂yₕ₁yₕ₂)⋅L)
                        GTensLConsDict[k₁,k₂,h₁,h₂] = [ρ[d₁*(k₁-1) + k₂, d₂*(h₁-1) + h₂]*x[transpose(LMB[i,:] + LMB[j,:] )] - x[transpose(LMB[i,:] + LMB[j,:] + eₖ₁ + eₕ₁ + eₖ₂ + eₕ₂)] for i in 1:nb_mon, j in 1:nb_mon ]
                    end
                end
            end
        end
    end
    return GTensLConsDict
end



function assembleSEPGTensLConsMat(d₁,d₂,GTensLConsDict)
    GTensLConMat  = []
    GTensLConCol = []
    for k₁ in 1:d₁
        for h₁ in 1:d₂
            for k₂ in 1:d₁
                for h₂ in 1:d₂
                    if k₂ == 1 && h₂ == 1
                        GTensLConCol = GTensLConsDict[k₁,k₂,h₁,h₂]
                    else
                        GTensLConCol = hcat(GTensLConCol, GTensLConsDict[k₁,k₂,h₁,h₂])
                    end
                end
            end
            println(size(GTensLConCol))
            if k₁ == 1 && h₁ == 1
                GTensLConMat  = GTensLConCol
            else
                GTensLConMat  = vcat(GTensLConMat,GTensLConCol)
            end
        end
    end
    return GTensLConMat
end

 k₁,k₂,h₁,h₂ = rand(1:d₁,(1,4))

GTensLConsDict = genSEPGTensLConsDict(ρ,t,x)
GTensLConMat = assembleSEPGTensLConsMat(d₁,d₂,GTensLConsDict)




"""M(G ⊗ L) ⪰ 0 constraints"""
function MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, isWeak=false)
    GρCon = []
    GρConRow = []


end
#@testset "genCPCons" begin
