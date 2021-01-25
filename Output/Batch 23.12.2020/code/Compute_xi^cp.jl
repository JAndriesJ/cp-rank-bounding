using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
using Dates # used in the output file naming

include("Make_xi^cp_constraints.jl")

"""This is where the ξₜᶜᵖ is calculated for matrix M """
function Computeξₜᶜᵖ(M,t = 1, isDag = false, GtensL = 0, isXX = false)
    n = size(M)[1]
    MomMatExp  = GenMomMatDict(n, t)
    MonBase = GenMon(n, t)
    nb_mon = size(MonBase)[1]

    ## Begin making the model
    model = Model(Mosek.Optimizer)
    # Define all variables that occur in the moment matrix.
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )

    # Build the moment matrix and constrain it to be PSD.
    MomMat = [x[transpose(MonBase[i,:] + MonBase[j,:])]   for i in 1:nb_mon, j in 1:nb_mon]
    @SDconstraint(model, MomMat >= zeros(nb_mon, nb_mon))

    # Localizing constraints (recall the form of g in the qudratic module)
    LMB = GenMon(n, t - 1)
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
        # for key in keys(preWeakGTensLConsDict)
        #     @SDconstraint(model, preWeakGTensLConsDict[key] >= zeros(LMB_nb, LMB_nb))
        # end

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


function  compξᶜᵖ(M,t)
    ξₜᶜᵖ           =  Computeξₜᶜᵖ(M, t, false, 0,false)
    ξₜᵩᶜᵖ          =  Computeξₜᶜᵖ(M, t, true,0,false)
    ξₜᵩweakTensᶜᵖ  =  Computeξₜᶜᵖ(M, t, true, 1,false)
    ξₜᵩTensᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, 2,false)
    ξₜᵩTensₓₓᶜᵖ    =  Computeξₜᶜᵖ(M, t, true, 2, true)
    return ξₜᶜᵖ, ξₜᵩᶜᵖ, ξₜᵩweakTensᶜᵖ,ξₜᵩTensᶜᵖ, ξₜᵩTensₓₓᶜᵖ
end










# function  Seach4BigCP(t)
#     for ind in 1:10
#         ξₜᶜᵖ        =  Computeξₜᶜᵖ(M, t, false, false, false)
#         ξₜᵩₔᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, true, false)
#     end
#     return  ξₜᶜᵖ, ξₜᵩₔᶜᵖ
# end

# "ξ₃ᶜᵖ,      ξ₃ᵩᶜᵖ,      ξ₃ᵩweak⊗ᶜᵖ,     ξ₃ᵩ⊗ᶜᵖ,    ξ₃ᵩ⊗ᶜᵖ + xᵢxⱼ "
