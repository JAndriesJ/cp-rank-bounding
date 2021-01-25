using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
using Dates # used in the output file naming
using Test # Utility module for testing if code works as intended.

include("xi^cp.jl")
include("Matrix_Examples_CP.jl") # Import matrices Examples we are interested in.


"""This is where the ξ_t is calculated for matrix M """
function Computeξₜᶜᵖ(M,t = 2, isDag = true, isGtensL = true, isXX = false)
    n = size(M)[1]
    MomMatExp  = GenMomMatDict(n, 2*t)
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
    LocMonBase = GenMon(n, t - 1)
    Loc_nb_mon  = size(LocMonBase)[1]

    LocConDict, DagConDict, GConDict, LocConXXDict = genCPCons(M,LocMonBase,x)

    ZOneByMon = zeros(Loc_nb_mon, Loc_nb_mon)
    ZMonByMon = zeros(Loc_nb_mon, Loc_nb_mon)

    for k in 1:n
        eₖ = LocMonBase[(k+1),:]
        #Second order Moment constraints
        @constraint(model, x[transpose(2*eₖ)]  == M[k,k])

        # Localizing g constriant
        @SDconstraint(model, LocConDict[(k,k)] >= ZMonByMon )
        # Localizing XX constriant
        if isXX
            @SDconstraint(model, LocConXXDict[(k,k)] >= ZMonByMon )
        end

        # Dagger constraints
        if isDag
            @constraint(model, DagConDict[(k,k)]  .>= ZOneByMon)
        end
        for h in (k+1):n
            eₕ = LocMonBase[(h+1),:]
            # Second order Moment constrain
            @constraint(model, x[transpose(eₖ + eₕ)]  == M[k,h])

            # Localizing g constriant
            @SDconstraint(model, LocConDict[(k,h)] >= ZMonByMon)
            # Localizing XX constriant
            if isXX
                @SDconstraint(model, LocConXXDict[(k,h)] >= ZMonByMon )
            end

            # Dagger constraints
            if isDag
                @constraint(model, DagConDict[(k,h)]  .>= ZOneByMon)
            end
        end
    end
    #G GenTensConst
    if isGtensL
        GCon = MakeG(LocConDict, GConDict, n)
        @SDconstraint(model, GCon >= zeros(Loc_nb_mon*n, Loc_nb_mon*n))
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

    return objective_value(model)
end

run_test = true
if run_test
    @testset "Computeξₜᶜᵖ" begin
        M  = M7 # your choices are: # M7 M7tilde M8tilde M9tilde M11tilde
        t  = 2
        ξ₂ᶜᵖ        =  Computeξₜᶜᵖ(M, t,false,false,false)
        ξ₂ᵩᶜᵖ       =  Computeξₜᶜᵖ(M, t, true, false, false)
        ξ₂ᵩₔᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, true, false)
        ξ₂ᵩₓₓᶜᵖ     =  Computeξₜᶜᵖ(M, t, true, true, true)
        @test ξ₂ᶜᵖ - 4.2182699647636985 == 0
        @test ξ₂ᵩᶜᵖ - 6.109713495002089 == 0
        @test ξ₂ᵩₔᶜᵖ - 6.956853350600587 == 0
        @test ξ₂ᵩₓₓᶜᵖ - 9.6386137347933 == 0
    end
end

"""Computes various ξₜ values at once"""
function  compξᶜᵖ(M,t)
    n²div4      =  floor(size(M)[1]^2 / 4)
    ξₜᶜᵖ        =  Computeξₜᶜᵖ(M, t, false, false, false)
    ξₜᵩᶜᵖ       =  Computeξₜᶜᵖ(M, t, true, false, false)
    ξₜᵩₔᶜᵖ      =  Computeξₜᶜᵖ(M, t, true, true, false)
    ξₜᵩₓₓᶜᵖ     =  Computeξₜᶜᵖ(M, t, true, true, true)
    return n²div4, ξₜᶜᵖ, ξₜᵩᶜᵖ, ξₜᵩₔᶜᵖ, ξₜᵩₓₓᶜᵖ
end

"""save computations to a .txt file """
function BatchCompξᶜᵖ(SavePath)
    #table1 =  Dict()
    open( SavePath * string(today()) * "RandomMats.txt", "w") do f
        write(f, "Matrix, " * "cp-rank, " * "n²/4, " * "ξ₂, " * "ξ₂ᵩ, " * "ξ₂ᵩₔ, "  *"ξ₂ᵩ + xᵢxⱼ, " *  "  \n")
        cp_ranks = [14 14 18 26 34]
        cp_matrices = [M7, M7tilde, M8tilde, M9tilde, M11tilde]
        cp_matrices_names = ["M7", "M7tilde", "M8tilde", "M9tilde", "M11tilde"]
        for ind in 1:5
            b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξᶜᵖ(cp_matrices[ind],2)
            write(f,  cp_matrices_names[ind]*","* string(cp_ranks[ind]) *",$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
        end

        counter = 1
        for  n in rand(6:9, (1,7))
            counter = counter + 1
            r = rand(4:11)
            M = genCPmatrix(n,r)
            Mtext = "$n × $n Rᶜᵖ_$r"
            b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξᶜᵖ(M,2)
            write(f,  Mtext*" , $r, $b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
            println("############################  $counter #############################  ")
        end

    end
end
BatchCompξᶜᵖ("C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Output\\")


# write(f,  "M7, 14,"* "$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
# b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξ(M7tilde,t)
# write(f,  "M7tilde, 14,"* "$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
# b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξ(M8tilde,t)
# write(f,  "M8tilde, 14,"* "$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
# b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξ(M9tilde,t)
# write(f,  "M9tilde, 14,"* "$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
# b, ξ₂, ξ₂ᵩ, ξ₂ᵩₔ, ξ₂ᵩₓₓ = compξ(M11tilde,t)
# write(f,  "M11tilde, 14,"* "$b, $ξ₂, $ξ₂ᵩ, $ξ₂ᵩₔ, $ξ₂ᵩₓₓ \n")
