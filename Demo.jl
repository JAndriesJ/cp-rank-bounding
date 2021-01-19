using Test # Utility module for testing if code works as intended.
#using JuMP # For the optimization frame work.
#using MosekTools # The solver that we use.

# Test Data generation:

include("C:\\Users\\andries\\all-my-codes\\matrix-menagerie\\mat_repo.jl")
import .mat_repo
include("moment_utils.jl")
include("constraints.jl")
include("compute.jl")


a,b,c,d,e = (0,0,0,0,0)
loadPath = "C:\\Users\\andries\\all-my-codes\\cp-rank-bounding\\Data\\CPmats\\M7.txt"
che = (n,t) -> binomial(n+1,t)
## matrices
if 1 == a
    @testset "loadMatfromtxt and saveMat2txt" begin
        M7 = mat_repo.loadMatfromtxt(loadPath)

        M = M7
        savePath = "C:\\Users\\andries\\all-my-codes\\cp-rank-bounding\\Data\\"
        saveName = "test-example"
        mat_repo.saveMat2txt(M, savePath, saveName )
    end

    @testset "testNN" begin
        M = rand(10,10)
        @test mat_repo.testNN(M) == true
    end

    @testset "testPSD" begin
        a = rand(10,1)
        M = a* tra(a)
        @test mat_repo.testPSD(M) == true
    end

    @testset "circmatrix" begin
        a = 1:4
        M = mat_repo.circmatrix(a)
        M_ref =   [4  3  2  1;
                    1  4  3  2;
                    2  1  4  3;
                    3  2  1  4]
        @test M == M_ref
    end

    @testset "makediagone" begin
    end

    @testset "genCPmatrix" begin
        A, a_ℓ= mat_repo.genRandCPmatrix(8,5)
        @test mat_repo.testNN(A)
        @test mat_repo.testPSD(A)
    end

    @testset "genSNNmatrix" begin
        A = mat_repo.genSNNmatrix(8)
        @test mat_repo.testNN(A)
        @test mat_repo.testPSD(A)
    end

    @testset "DominateDiagonal(A,λ)" begin
    end

    @testset "genDDSNNmatrix" begin
    end
end

## moments
# This code checks if the variables representing the exponents of moments are coded as intended
if 1 == b

    # Moments
    @testset "make_mon_expo(n::Int,t::Int, isLeq::Bool = true)" begin
        for n ∈ 3:9, t ∈ 0:4
            @test size(make_mon_expo_arr(n,t))[2] == n
            MonBase = make_mon_expo(n,t,true)
            #@test length(MonBase) == che(n,t)
        end
        for n ∈ 3:9, t ∈ 0:3
            MonBase = make_mon_expo(n,t,true)
            @test length(MonBase) == binomial(n +t,t)
        end

        MonBase = make_mon_expo(2,3,false)
        @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]

        MonBase = make_mon_expo(2,4,false)
        @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]
    end

    @testset "get_std_base_vec(n::Int,k::Int)" begin
        n =7
        for k in 1:n
            e = get_std_base_vec(n,k)
            @test length(e) == n
            @test maximum(e) == 1
            @test minimum(e) == 0
            @test sum(e) == 1
        end
     end

    @testset "get_mon_index(B,α)" begin
        MonBase =  make_mon_expo(2,4,false)
        α = [2  2]
        @test  get_mon_index(MonBase, α) == 3

        MonBase = make_mon_expo(2,3)
        α = [2.0 1.0]
        @test  get_mon_index(MonBase, α) == 8
    end

    @testset "make_mom__expo_mat_dict(n::Int,t::Int)" begin
        n = 2
        t = 2
        MonExp = make_mom_expo_mat_dict(n,t)
        MonBase = make_mon_expo(n,t)
        for key in keys(MonExp)
            for val in MonExp[key]
                @test   key == MonBase[val[1]] + MonBase[val[2]]
            end
        end
    end
end

# This code checks if the constraints are implemented as intended.
if 1 == c
    A = mat_repo.loadMatfromtxt(loadPath)
    t = 2
    n = size(A)[1]
    Lx = make_dummy_var(n,t)

    @testset "make_loc_con" begin
        loc_con = make_loc_con(A,t,Lx)
        @test length(keys(loc_con)) == che(n,t)
    end

    @testset "make_dag_con" begin
        dag_con = make_dag_con(A,t,Lx)
        @test length(keys(dag_con)) == n*(n+1)/2 + 1
    end

    @testset "make_xx_con" begin
        xx_con = make_xx_con(A,t,Lx)
        @test length(keys(xx_con)) == n*(n+1)/4 + n ### ???
    end

    @testset "make_weakG_con" begin
        weakG_con = make_weakG_con(A,t,Lx)
        length(keys(weakG_con)) == t
        for key in keys(weakG_con)
            size(weakG_con[key]) == (n^key,n^key)
        end
    end

    @testset "make_G_con" begin
        #LocConDict         = genCP_localizing_Constraints(A,MonBaseₜ₋₁,Lx)
        #GTensLConsDict     = MakeGTensLConsMat1(LocConDict, A, MonBaseₜ₋₁,Lx)
        G_con     = make_G_con(A,t,Lx)
        size(G_con) == (n + n^2,n + n^2)
    end
end



# Test Computeξₜᶜᵖ for example matrices.
if 1 == d
    @testset "Computeξₜᶜᵖ" begin
        cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
        loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[3]
        A = mat_repo.loadMatfromtxt(loadPath)

        t  = 2
        ξ₂ᶜᵖ           = Computeξₜᶜᵖ(A, t, false,0,false)
        @test ξ₂ᶜᵖ           == 4.2183

        ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(A, t, true, 0,false)
        @test ξ₂ᵩᶜᵖ          == 6.2388

        ξ₂ₓₓᶜᵖ         = Computeξₜᶜᵖ(A, t, false, 0, true)
        @test ξ₂ₓₓᶜᵖ         == 5.0581

        ξ₂weakTensᶜᵖ  = Computeξₜᶜᵖ(A, t, false, 1,false)
        @test ξ₂weakTensᶜᵖ   == 6.907

        ξ₂Tensᶜᵖ      = Computeξₜᶜᵖ(A, t, false, 2,false)
        @test ξ₂Tensᶜᵖ       == 6.8033

        ξ₂ᵩTensᶜᵖ      = Computeξₜᶜᵖ(A, t, true, 2,false)
        @test ξ₂ᵩTensᶜᵖ      == 6.9999

        ξₜᵩweakTensᶜᵖ  = Computeξₜᶜᵖ(A, t, true, 1,false)
        @test ξ₂ᵩweakTensᶜᵖ  == 7.0

        ξ₂ᵩTensₓₓᶜᵖ    = Computeξₜᶜᵖ(A, t, true, 2, true)
        @test ξ₂ᵩTensₓₓᶜᵖ    == 9.8757


    end
end


## moments

## load a matrix:

cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[5]
A = mat_repo.loadMatfromtxt(loadPath)

t  = 2
ξ₂ᶜᵖ           = Computeξₜᶜᵖ(A, t, false,0,false)
#@test ξ₂ᶜᵖ           == 4.2183

ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(A, t, true, 0,false)
#@test ξ₂ᵩᶜᵖ          == 6.2388

ξ₂ₓₓᶜᵖ         = Computeξₜᶜᵖ(A, t, false, 0, true)
#@test ξ₂ₓₓᶜᵖ         == 5.0581



ξ₂weakGᶜᵖ  = Computeξₜᶜᵖ(A, t, false, 1,false)
#@test ξ₂weakTensᶜᵖ   == 6.907

ξ₂Gᶜᵖ      = Computeξₜᶜᵖ(A, t, false, 2,false)
#@test ξ₂Tensᶜᵖ       == 6.8033


##  https://jump.dev/JuMP.jl/dev/reference/solutions/


# cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
# loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[3]
# M = mat_repo.loadMatfromtxt(loadPath)
# t = 2
#
#
# n = size(M)[1]
# MomMatExp  = make_mom__expo_mat_dict(n, t)
# MonBase = make_mon_expo(n, t)
# nb_mon = size(MonBase)[1]
# model = Model(Mosek.Optimizer)
#
# list_of_keys = [key for key in keys(MomMatExp) ]
# @variable(model, Lx[list_of_keys] )
#
#
# ar  = genCPweakGTensLCons(M,t,Lx)
# nar = MakeGTensLConsMat(M,t,Lx)
# kar = MakeGTensLConsMat2(M,t,Lx)
#
#
# for i in 1:56
#     for j in 1:56
#         if (nar[i,j] != kar[i,j])
#             println("$i,$j")
#         end
#     end
# end
#
#
#
# LMBexp_00         = make_mon_expo_mat(n,0,false) #exponents of [x]₌₀[x]₌₀ᵀ
# LMBexp_10         = make_mon_expo_mat(n,(1,0),false) #exponents of [x]₌₁[x]₌₀ᵀ
# LMBexp_01         = make_mon_expo_mat(n,(0,1),false)
# LMBexp_11         = make_mon_expo_mat(n,1,false) #exponents of [x]₌₁[x]₌₁ᵀ
#
#
# LMBexp_leq11      = vcat(hcat(LMBexp_00,LMBexp_01),hcat(LMBexp_10,LMBexp_11 ))
