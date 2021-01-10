cl = clearconsole
using Test # Utility module for testing if code works as intended.
#using JuMP # For the optimization frame work.
#using MosekTools # The solver that we use.

# Test Data generation:

include("C:\\Users\\andries\\all-my-codes\\matrix-menagerie\\mat_repo.jl")
import .mat_repo

# include("MatUtils.jl")
# include("Generate_data.jl")

include("moment_utils.jl")
include("Make_xi^cp_constraints.jl")
include("Compute_xi^cp.jl")


a,b,c,d,e = (0,0,0,0,1)
loadPath = "C:\\Users\\andries\\all-my-codes\\cp-rank-bounding\\Data\\CPmats\\M7.txt"

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
        M = a* tr(a)
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
# This code checks if the variables representing the moments are coded as intended

if 1 == b
    # Moments

    @testset "GenMon" begin
        MonBase = GenMon(2,3,false)
        @test MonBase ==  [3.0  0.0;
                           2.0  1.0;
                           1.0  2.0;
                           0.0  3.0]

        MonBase = GenMon(2,4,false)
        @test MonBase ==  [4  0;
                           3  1;
                           2  2;
                           1  3;
                           0  4]
    end

    @testset "standardBase(n,k)" begin
        n =7
        for k in 1:n
            e = standardBase(n,k)
            @test length(e) == n
            @test maximum(e) == 1
            @test minimum(e) == 0
            @test sum(e) == 1
        end
     end



    @testset "GetMonIndex" begin
        MonBase =  GenMon(2,4,false)
        α = [2  2]
        @test  GetMonIndex(MonBase, α) == 3

        MonBase = GenMon(2,3)
        α = [2.0 1.0]
        @test  GetMonIndex(MonBase, α) == 8
    end

    @testset "makeArray" begin
       #
    end

    @testset "GenMomMatDict" begin
        n = 2
        t = 2
        MonExp = GenMomMatDict(n,t)
        MonBase = GenMon(n,t)
        for key in keys(MonExp)
            for val in MonExp[key]
                @test   makeArray(key) == makeArray(MonBase[val[1],:] + MonBase[val[2],:])
            end
        end
    end
end

# This code checks if the constraints are implemented as intended.

if 1 == c
    M = mat_repo.loadMatfromtxt(loadPath)
    t = 2
    n = size(M)[1]
    MomMatExp  = GenMomMatDict(n, t)
    MonBase = GenMon(n, t)
    model = Model(Mosek.Optimizer)
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )
    LMB = GenMon(n, t-1)


    @testset "genCP_localizing_Constriaints" begin
        LocConDict = genCP_localizing_Constriaints(M,LMB,x)
    end

    @testset "genCP_dagger_Constriaints" begin
        DagConDict = genCP_dagger_Constriaints(M,LMB,x)
        @test length(keys(DagConDict)) == n*(n+1)/2
    end

    @testset "genCP_XX_Constriaints" begin
        XXConDict = genCP_XX_Constriaints(M,LMB,x)
    end

    @testset "genCPpreWeakGTensLCons" begin
        LocConDict         = genCP_localizing_Constriaints(M,LMB,x)
        WeakGTensLConsDict = genCPweakGTensLCons(M,LMB,x)
        GTensLConsDict     = MakeGTensLConsMat(LocConDict, WeakGTensLConsDict, n)
    end
end

# Test Computeξₜᶜᵖ for example matrices.
if 1 == d
    @testset "Computeξₜᶜᵖ" begin
        M  = mat_repo.loadMatfromtxt(loadPath)
        t  = 2
        ξ₂ᶜᵖ           = Computeξₜᶜᵖ(M, t, false,0,false)
        ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(M, t, true, 0,false)
        ξₜᵩweakTensᶜᵖ  = Computeξₜᶜᵖ(M, t, true, 1,false)
        ξ₂ᵩTensᶜᵖ      = Computeξₜᶜᵖ(M, t, true, 2,false)
        ξ₂ᵩTensₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, true, 2, true)
        @test ξ₂ᶜᵖ          - 4.2183 == 0
        @test ξ₂ᵩᶜᵖ         - 6.1097 == 0
        @test ξₜᵩweakTensᶜᵖ - 6.9569 == 0
        @test ξ₂ᵩTensᶜᵖ     - 6.9569 == 0
        @test ξ₂ᵩTensₓₓᶜᵖ   - 9.6389 == 0
    end
end




## moments

## load a matrix:




#if 1 == 1

cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
# rand_cp_mats = ["1 R n7 r6.txt"    "14 R n6 r4.txt"   "19 R n6 r15.txt"  "5 R n8 r26.txt"]
#                 # '10 R n9 r28.txt'  '15 R n8 r7.txt'   '2 R n6 r23.txt'   '6 R n7 r6.txt'
#                 # '11 R n7 r21.txt'  '16 R n6 r16.txt'  '20 R n6 r22.txt'  '7 R n7 r23.txt'
#                 # '12 R n9 r28.txt'  '17 R n9 r24.txt'  '3 R n6 r29.txt'   '8 R n8 r11.txt'
#                 # '13 R n6 r19.txt'  '18 R n9 r17.txt'  '4 R n8 r30.txt'   '9 R n9 r19.txt']
#
#
loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[3]
M = mat_repo.loadMatfromtxt(loadPath)
t  = 2

# ξ₂ᶜᵖ           = Computeξₜᶜᵖ(M, t, false,0,false)
# ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(M, t, true, 0,false)

n = size(M)[1]
MomMatExp  = GenMomMatDict(n, t)
MonBase = GenMon(n, t)
model = Model(Mosek.Optimizer)
list_of_keys = [key for key in keys(MomMatExp) ]
@variable(model, x[list_of_keys] )

t_temp = 1
mon_eq_mat_ex  = gen_mon_eq_mat(n,t_temp)
mon_eq_mat     = index_to_var(x, mon_eq_mat_ex)
K = size(mon_eq_mat)[1]
g_mon_eq_mat = Dict()
i = 1
eᵢ = standardBase(n,i)
j = 1
eⱼ = standardBase(n,j)
B = repeat( [eᵢ + eⱼ] , inner = (1,1), outer = (K,K))
# g_mon_eq_mat[(i,j)] = M[i,j]*mon_eq_mat - index_to_var(x, mon_eq_mat_ex + B)

g_mon_eq_mat   = gen_g_mon_eq_mat(M,mon_eq_mat,mon_eq_mat_ex,x)
weakGTensLCons[t_temp] = g_mon_eq_mat


WeakGTensLConsDict = genCPweakGTensLCons(M,t,x)

ξ₂weakTensᶜᵖ   = Computeξₜᶜᵖ(M, t, false, 1,false)
ξ₂ᵩweakTensᶜᵖ  = Computeξₜᶜᵖ(M, t, true, 1,false)
#
# # #
# println("ξ₂ᶜᵖ: $ξ₂ᶜᵖ ,ξ₂ᵩᶜᵖ:  $ξ₂ᵩᶜᵖ,ξ₂weakTensᶜᵖ: $ξ₂weakTensᶜᵖ ,ξ₂ᵩweakTensᶜᵖ: $ξ₂ᵩweakTensᶜᵖ")
    # ξ₂ᵩTensᶜᵖ      = Computeξₜᶜᵖ(M, t, true, 2,false)
    # ξ₂Tensᶜᵖ       = Computeξₜᶜᵖ(M, t, false, 2,false)
    # ξ₂ᵩTensₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, true, 2, true)
    # ξ₂ᵩₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, true, 0, true)
    # ξ₂ₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, false, 0, true)





    # n = size(M)[1]
    # MomMatExp  = GenMomMatDict(n, t)
    # MonBase = GenMon(n, t)
    # model = Model(Mosek.Optimizer)
    # list_of_keys = [key for key in keys(MomMatExp) ]
    # @variable(model, x[list_of_keys] )



    #DagConDict = genCP_dagger_Constriaints(M,t,x)



    # nar =  [[0  0  0  0  0  0  0], [1  0  0  0  0  0  0] , [0  1  0  0  0  0  0] , [0  0  1  0  0  0  0]]
    # nar = reshape(nar ,2,2)
    # index_to_var(x, nar)
    #
    #
    # poes = reshape(tr( [LMB[i,:] + LMB[j,:]  for i in 1:8 for j in 1:8 ]), 8,8)
    #
    # nar = [[0 0] ,[1 0], [0 1], [1 1]]
    # nar = reshape(nar,2,2)
    # nar = nar + repeat([[1 0]], inner=(1, 1), outer=(2, 2))
#end

# fielName_lst = ["(1)-6 Rcp_6","(11)-7 Rcp_6","(21)-8 Rcp_21", "(32)-9 Rcp_5"]
# loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPrandMat\\"*fielName_lst[3]*".txt"
# M = loadMatfromtxt(loadPath)
# t = 2
# n = size(M)[1]
# MomMatExp  = GenMomMatDict(n, t)
# MonBase = GenMon(n, t)
# model = Model(Mosek.Optimizer)
# list_of_keys = [key for key in keys(MomMatExp) ]
# @variable(model, x[list_of_keys] )
# LMB = GenMon(n, t-1)
#
#
# LocConDict            = genCP_localizing_Constriaints(M,LMB,x)
# DagConDict            = genCP_dagger_Constriaints(M,LMB,x)
# LocConXXDict          = genCP_XX_Constriaints(M,LMB,x)
#
# weakGTensLConsMat     = genCPweakGTensLCons(M,LMB,x)
# GTensLConsMat         = MakeGTensLConsMat(LocConDict, weakGTensLConsMat, n)
#
#
#
# ξₜᵩᶜᵖ         =  Computeξₜᶜᵖ(M, t, true, 0,false)
# ξₜweakTensᶜᵖ  = Computeξₜᶜᵖ(M, t, false, 1,false)
#
# println("ξₜᵩᶜᵖ: $ξₜᵩᶜᵖ     ξₜweakTensᶜᵖ: $ξₜweakTensᶜᵖ")
