using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
using Dates # used in the output file naming

include("Make_xi^cp_constraints.jl")

"""This is where the ξₜᶜᵖ is calculated for matrix M """
function Computeξₜᶜᵖ(M,t , isDag , GtensL, isXX )
    n = size(M)[1]
    MomMatExp  = make_mom__expo_mat_dict(n, t)
    MonBase = make_mon_expo(n, t)
    nb_mon = size(MonBase)[1]

    ## Begin making the model
    model = Model(Mosek.Optimizer)
    # Define all variables that occur in the moment matrix.
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )

## Build the moment matrix and constrain it to be PSD.
    MomMat = [x[transpose(MonBase[i,:] + MonBase[j,:])]   for i in 1:nb_mon, j in 1:nb_mon]
    @SDconstraint(model, MomMat >= zeros(nb_mon, nb_mon))

## Localizing constraints (recall the form of g in the qudratic module)
    LMB = make_mon_expo(n, t - 1)
    LMB_nb  = size(LMB)[1]


    ZOneByMon = zeros(1, LMB_nb)
    ZMonByMon = zeros(LMB_nb, LMB_nb)


## Second order Moment constraints
    for k in 1:n
        eₖ = get_standard_base(n,k)
        @constraint(model, x[transpose(2*eₖ)]  == M[k,k])
        for h in (k+1):n
            eₕ  = get_standard_base(n,h)
            @constraint(model, x[transpose(eₖ + eₕ)]  == M[k,h])
        end
    end

## Localizing g constraint
    LocConDict = genCP_localizing_Constraints(M,LMB,x)
    for k in 1:n
        eₖ = get_standard_base(n,k)
        @SDconstraint(model, LocConDict[(k,k)] >= ZMonByMon )
        for h in (k+1):n
            eₕ  = get_standard_base(n,h)
            @SDconstraint(model, LocConDict[(k,h)] >= ZMonByMon)
        end
    end

## Dagger constraints
    if isDag
        println("----------------Dagger constraints are active")
        DagConDict  = genCP_dagger_Constraints(M,t,x)
                # MISTAKE
        # for k in 1:n
        #     eₖ = standardBase(n,k)
        #     for i in 1:LMB_nb
        #         @constraint(model, DagConDict[(k,k)][i]  >= 0)
        #     end
        #     for h in (k+1):n
        #         eₕ  = standardBase(n,h)
        #         for i in 1:LMB_nb
        #             @constraint(model, DagConDict[(k,h)][i]  >= 0)
        #         end
        #     end
        # end

        for k in 1:n
            eₖ = get_standard_base(n,k)
            @constraint(model, DagConDict[(k,k)]  .>= 0)
            for h in (k+1):n
                eₕ  = get_standard_base(n,h)
                @constraint(model, DagConDict[(k,h)]  .>= 0)
            end
        end
    end

    # Localizing XX constraint
    if isXX
        println("----------------XX constraints are active")
        for k in 1:n
            eₖ = get_standard_base(n,k)
                LocConXXDict          = genCP_XX_Constraints(M,LMB,x)
                @SDconstraint(model, LocConXXDict[(k,k)] >= ZMonByMon )
            for h in (k+1):n
                eₕ = get_standard_base(n,h)
                @SDconstraint(model, LocConXXDict[(k,h)] >= ZMonByMon )
            end
        end
    end

    #G GenWeakTensConst
    if GtensL == 0

    elseif GtensL == 1
        println("----------------Weak tensor-constraints are active")
        WeakGTensLConsDict = genCPweakGTensLCons(M,t,x)
        @SDconstraint(model, WeakGTensLConsDict[0] >= zeros(size(WeakGTensLConsDict[0])))
        @SDconstraint(model, WeakGTensLConsDict[1] >= zeros(size(WeakGTensLConsDict[1])))
        # for key in keys(WeakGTensLConsDict)
           # @SDconstraint(model, WeakGTensLConsDict[key] >= zeros(size(WeakGTensLConsDict[key])))
        # end

    elseif GtensL == 2
        println("----------------Tensor-constraints are active")
        GTensLConsMat       = MakeGTensLConsMat1(LocConDict, M, LMB,x)
        #GTensLConsMat       = MakeGTensLConsMat(M,t,x)
        # for i in 1:8
        #     for j in 1:8
        #             if GTensLConsMat[i,j] != GTensLConsMat[i,j]
        #                 println("$i,$j")
        #             end
        #     end
        # end

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
    if string(primal_status(model)) == "NO_SOLUTION"
        return NaN
    else
        return round(objective_value(model), digits = 4)
    end

end


# function  compξᶜᵖ(M,t)
#     ξₜᶜᵖ           =  Computeξₜᶜᵖ(M, t, false, 0,false)
#     ξₜᵩᶜᵖ          =  Computeξₜᶜᵖ(M, t, true,0,false)
#     ξₜᵩweakTensᶜᵖ  =  Computeξₜᶜᵖ(M, t, true, 1,false)
#     ξₜweakTensᶜᵖ  =  Computeξₜᶜᵖ(M, t, false, 1,false)
#     ξₜᵩTensᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, 2,false)
#     ξₜTensᶜᵖ      =  Computeξₜᶜᵖ(M, t, false, 2,false)
#     ξₜᵩTensₓₓᶜᵖ    =  Computeξₜᶜᵖ(M, t, true, 2, true)
#     return ξₜᶜᵖ, ξₜᵩᶜᵖ, ξₜweakTensᶜᵖ,ξₜᵩweakTensᶜᵖ,ξₜTensᶜᵖ,ξₜᵩTensᶜᵖ, ξₜᵩTensₓₓᶜᵖ
# end
#
# function  compξᶜᵖ_1(M,t)
#     ξₜᶜᵖ           =  Computeξₜᶜᵖ(M, t, false, 0,false)
#     ξₜᵩᶜᵖ          =  Computeξₜᶜᵖ(M, t, true,0,false)
#     ξₜweakTensᶜᵖ  =  Computeξₜᶜᵖ(M, t, false, 1,false)
#     ξₜTensᶜᵖ      =  Computeξₜᶜᵖ(M, t, false, 2,false)
#     return ξₜᶜᵖ, ξₜᵩᶜᵖ, ξₜweakTensᶜᵖ,ξₜTensᶜᵖ
# end

 # ξₜᵩweakTensᶜᵖ  =  Computeξₜᶜᵖ(M, t, true, 1,false)
 # ξₜᵩTensᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, 2,false)
 # ξₜᵩTensₓₓᶜᵖ    =  Computeξₜᶜᵖ(M, t, true, 2, true)





# function  Seach4BigCP(t)
#     for ind in 1:10
#         ξₜᶜᵖ        =  Computeξₜᶜᵖ(M, t, false, false, false)
#         ξₜᵩₔᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, true, false)
#     end
#     return  ξₜᶜᵖ, ξₜᵩₔᶜᵖ
# end

# "ξ₃ᶜᵖ,      ξ₃ᵩᶜᵖ,      ξ₃ᵩweak⊗ᶜᵖ,     ξ₃ᵩ⊗ᶜᵖ,    ξ₃ᵩ⊗ᶜᵖ + xᵢxⱼ "
