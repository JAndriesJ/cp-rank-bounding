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


a,b,c,d,e = (0,0,0,0,0)
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

# This code checks if the variables representing the exponents of moments are coded as intended
if 1 == b
    # Moments
    @testset "GenMon" begin
        MonBase = make_mon_expo(2,3,false)
        @test MonBase ==  [3.0  0.0;
                           2.0  1.0;
                           1.0  2.0;
                           0.0  3.0]

        MonBase = make_mon_expo(2,4,false)
        @test MonBase ==  [4  0;
                           3  1;
                           2  2;
                           1  3;
                           0  4]
    end

    @testset "standardBase(n,k)" begin
        n =7
        for k in 1:n
            e = get_standard_base(n,k)
            @test length(e) == n
            @test maximum(e) == 1
            @test minimum(e) == 0
            @test sum(e) == 1
        end
     end

    @testset "GetMonIndex" begin
        MonBase =  make_mon_expo(2,4,false)
        α = [2  2]
        @test  get_mon_index(MonBase, α) == 3

        MonBase = make_mon_expo(2,3)
        α = [2.0 1.0]
        @test  get_mon_index(MonBase, α) == 8
    end

    @testset "makeArray" begin
       #
    end

    @testset "GenMomMatDict" begin
        n = 2
        t = 2
        MonExp = make_mom__expo_mat_dict(n,t)
        MonBase = make_mon_expo(n,t)
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

    MomMatExp  = make_mom__expo_mat_dict(n, t)
    MonBaseₜ    = make_mon_expo(n, t)
    model      = Model(Mosek.Optimizer)
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model,Lx[list_of_keys] )

    MonBaseₜ₋₁ = make_mon_expo(n, t-1)


    @testset "genCP_localizing_Constriaints" begin
        LocConDict = genCP_localizing_Constriaints(M,MonBaseₜ₋₁,Lx)
    end

    @testset "genCP_dagger_Constriaints" begin
        DagConDict = genCP_dagger_Constriaints(M,t,Lx)
        @test length(keys(DagConDict)) == n*(n+1)/2 + 1
    end

    @testset "genCP_XX_Constriaints" begin
        XXConDict = genCP_XX_Constriaints(M,MonBaseₜ₋₁,Lx)
    end

    @testset "genCPweakGTensLCons" begin
        WeakGTensLConsDict = genCPweakGTensLCons(M,t,Lx)
    end

    @testset "genCPGTensLCons" begin
        LocConDict = genCP_localizing_Constriaints(M,MonBaseₜ₋₁,Lx)
        GTensLConsDict     = MakeGTensLConsMat(LocConDict, M, MonBaseₜ₋₁,Lx)
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
        ξ₂ₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, false, 0, true)
        # @test ξ₂ᶜᵖ          - 4.2183 == 0
        # @test ξ₂ᵩᶜᵖ         - 6.1097 == 0
        # @test ξₜᵩweakTensᶜᵖ - 6.9569 == 0
        # @test ξ₂ᵩTensᶜᵖ     - 6.9569 == 0
        # @test ξ₂ᵩTensₓₓᶜᵖ   - 9.6389 == 0
    end
end








## moments

## load a matrix:
cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[3]
M = mat_repo.loadMatfromtxt(loadPath)
t = 2
n = size(M)[1]
# include("Compute_xi^cp.jl")


MomMatExp    = make_mom__expo_mat_dict(n, t)
MonBaseₜ     = make_mon_expo(n, t)
model        = Model(Mosek.Optimizer)
list_of_keys = [key for key in keys(MomMatExp) ]
@variable(model,Lx[list_of_keys] )


A             = MakeGTensLConsMat(M,t,Lx)

LMB           = make_mon_expo(n, t - 1)
LocConDict    = genCP_localizing_Constriaints(M,LMB,Lx)
B             = MakeGTensLConsMat1(LocConDict, M, LMB,Lx)
for i in 1:8
    for j in 1:8
            if A[i,j] != B[i,j]
                println("$i,$j")
            end
    end
end
for j in 1:8
    if A[:,j] != B[:,j]
        println(":,$j")
    end
end


A == B




# ξₜweakTensᶜᵖ  = Computeξₜᶜᵖ(M, t, false, 1,false)
# ξₜTensᶜᵖ  = Computeξₜᶜᵖ(M, t, false, 2,false)




#if 1 == 1

# M = mat_repo.loadMatfromtxt(loadPath)
# t  = 2
# n = size(M)[1]
#

#
# MonBaseₜ₋₁ = make_mon_expo(n, t-1)
#
#
#
# LocConDict         = genCP_localizing_Constriaints(M,MonBaseₜ₋₁,Lx)
# DagConDict         = genCP_dagger_Constriaints(M,t,Lx)
# XXConDict          = genCP_XX_Constriaints(M,MonBaseₜ₋₁,Lx)
#
#
# GTensLConsDict     = MakeGTensLConsMat(LocConDict, M, MonBaseₜ₋₁)
