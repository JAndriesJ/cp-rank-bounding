using Test # Utility module for testing if code works as intended.

# Data generation:
# MatUtils.jl and Generate_data.jl
include("MatUtils.jl")
include("Generate_data.jl")
testMatGen = false
if testMatGen
    @testset "saveMat2txt" begin
        M = M7
        savePath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\"
        saveName = "Example"
        saveMat2txt(M, savePath, saveName )
    end

    @testset "loadMatfromtxt" begin
        loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPrandMat\\6 Rcp_14-11-56-04.txt"
        M = loadMatfromtxt(loadPath)
    end

    @testset "testNN" begin
        M = rand(10,10)
        @test testNN(M) == true
    end

    @testset "testPSD" begin
        a = rand(10,1)
        M = a* transpose(a)
        @test testPSD(M) == true
    end


    @testset "circmatrix" begin
        a = 1:4
        M = circmatrix(a)
        M_ref =   [4  3  2  1;
                    1  4  3  2;
                    2  1  4  3;
                    3  2  1  4]
        @test M == M_ref
    end

    @testset "makediagone" begin

    end

    @testset "genCPmatrix" begin
        A, a_ℓ= genCPmatrix(8,5)
        @test testNN(A)
        @test testPSD(A)
    end


    @testset "genSNNmatrix" begin
        A = genSNNmatrix(8)
        @test testNN(A)
        @test testPSD(A)
    end


    @testset "DominateDiagonal(A,λ)" begin

    end


    @testset "genDDSNNmatrix" begin

    end
end

testMoments = false
if testMoments
    # Moments
    include("moment_utils.jl")

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

    # @testset "GenMomMatDictComplex" begin
    #     n = 2
    #     t = 2
    #     MonExp = GenMomMatDict(n,t*2)
    #     MonBase = GenMon(n,t)
    #     for key in keys(MonExp)
    #         for val in MonExp[key]
    #             @test   makeArray(key) == makeArray(MonBase[val[1],:] + MonBase[val[2],:])
    #         end
    #     end
    # end
end



# xi^cp
using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
include("Compute_xi^cp.jl")


if false

    M = M7
    t = 2
    n = size(M7)[1]
    MomMatExp  = GenMomMatDict(n, t)
    MonBase = GenMon(n, t)
    model = Model(Mosek.Optimizer)
    list_of_keys = [key for key in keys(MomMatExp) ]
    @variable(model, x[list_of_keys] )
    LMB = GenMon(n, 1)


    @testset "genCP_localizing_Constriaints" begin
        LocConDict = genCP_localizing_Constriaints(M,LMB,x)
    end

    @testset "genCP_dagger_Constriaints" begin
        DagConDict = genCP_dagger_Constriaints(M,LMB,x)
    end

    @testset "genCP_XX_Constriaints" begin
        XXConDict = genCP_XX_Constriaints(M,LMB,x)
    end

    @testset "genCPpreWeakGTensLCons" begin
        preWeakGTensLConsDict = genCPpreWeakGTensLCons(M,LMB,x)
    end


    @testset "genCPpreWeakGTensLCons" begin
        LocConDict = genCP_localizing_Constriaints(M,LMB,x)
        preWeakGTensLConsDict = genCPpreWeakGTensLCons(M,LMB,x)
        weakGTensLConsDict = MakeGTensLConsMat(LocConDict, preWeakGTensLConsDict, n, false)
    end

end

# Computeξₜᶜᵖ
if false
    @testset "Computeξₜᶜᵖ" begin
        M  = M7 # your choices are: # M7 M7tilde M8tilde M9tilde M11tilde
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
