using JuMP # For the optimization frame work.
using LinearAlgebra # For the matrix manipulations and preprocessing.
using Random # for the generation of a random cp Matrix
using Test # Utility module for testing if code works as intended.

include("xi^SEP.jl")
include("Matrix_Examples_SEP.jl") # Import matrices Examples we are interested in.

"""This is where the ξₜᶜᵖ is calculated for matrix M """
function Computeξₜᶜᵖ(ρ,t = 1, GtensL = 0)
    n = Int(2*sqrt(size(ρ)[1]))

    MomMatExp  = GenMomMatDict(n, t)
    MonBase = GenMon(n, t)
    nb_mon = size(MonBase)[1]
    model = Model(Mosek.Optimizer)
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )

    # Build the moment matrix and constrain it to be PSD.
    MomMat = [x[transpose(MonBase[i,:] + MonBase[j,:])]   for i in 1:nb_mon, j in 1:nb_mon]
    @SDconstraint(model, MomMat >= zeros(nb_mon, nb_mon))

    # Localizing constraints (recall the form of g in the qudratic module)

    LMB_nb  = size(LMB)[1]


    ZOneByMon = zeros(LMB_nb, LMB_nb)
    ZMonByMon = zeros(LMB_nb, LMB_nb)

    #Second order Moment constraints
    for k in 1:n
        eₖ = standardBase(n,k)
        @constraint(model, x[transpose(2*eₖ)]  == M[k,k])
        for h in (k+1):n
            eₕ  = standardBase(n,h)
            @constraint(model, x[transpose(eₖ + eₕ)]  == M[k,h])
        end
    end



    

    # Localizing g constriant
    LocConDict = genCP_localizing_Constriaints(M,LMB,x)
    for k in 1:n
        eₖ = standardBase(n,k)
        @SDconstraint(model, LocConDict[(k,k)] >= ZMonByMon )
        for h in (k+1):n
            eₕ  = standardBase(n,h)
            @SDconstraint(model, LocConDict[(k,h)] >= ZMonByMon)
        end
    end

    # Dagger constraints
    if isDag
        DagConDict            = genCP_dagger_Constriaints(M,LMB,x)
        for k in 1:n
            eₖ = standardBase(n,k)
            @constraint(model, DagConDict[(k,k)]  .>= ZOneByMon)
            for h in (k+1):n
                eₕ  = standardBase(n,h)
                @constraint(model, DagConDict[(k,h)]  .>= ZOneByMon)
            end
        end
    end

    # Localizing XX constriant
    if isXX
        for k in 1:n
            eₖ = standardBase(n,k)

                LocConXXDict          = genCP_XX_Constriaints(M,LMB,x)
                @SDconstraint(model, LocConXXDict[(k,k)] >= ZMonByMon )
            for h in (k+1):n
                eₕ = standardBase(n,h)
                @SDconstraint(model, LocConXXDict[(k,h)] >= ZMonByMon )
            end
        end
    end

    #G GenWeakTensConst
    if GtensL == 0

    elseif GtensL == 1
        preWeakGTensLConsDict = genCPpreWeakGTensLCons(M,LMB,x)
        weakGTensLConsMat     = MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, true)
        @SDconstraint(model, weakGTensLConsMat >= zeros(LMB_nb*n, LMB_nb*n))
    elseif GtensL == 2
        preWeakGTensLConsDict = genCPpreWeakGTensLCons(M,LMB,x)
        GTensLConsMat       = MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, false)
        @SDconstraint(model, GTensLConsMat >= zeros(LMB_nb*n, LMB_nb*n))
    end

    #  Set objective
    @objective(model, Min, x[transpose(zeros(n))] )
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    # println("Variables: ", "x = ",value(x))

    return round(objective_value(model), digits = 4)
end
