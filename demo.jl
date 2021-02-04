using Test # Utility module for testing if code works as intended.
using BenchmarkTools
#using JuMP # For the optimization frame work.
#using MosekTools # The solver that we use.

    include("mat_IO.jl")
    include("moment_utils.jl")
    include("constraints.jl")
    include("compute.jl")
    script_dir = @__DIR__


function quick_load_mat(n::Integer)
    cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
    loadPath = script_dir*"\\Data\\CPmats\\"*cp_mats[n]
    return loadMatfromtxt(loadPath)
end

A = quick_load_mat(6)
# A = rand(4,4)
# A = A*tra(A)

n = size(A)[1]
t  = 2

Lx_1,model_ξ₂_1ᶜᵖ  = Computeξₜᶜᵖ(A, t, false,1,false)
objective_value(model_ξ₂_1ᶜᵖ)

@show relative_gap(model_ξ₂_1ᶜᵖ)

function run_tests()
    a,b,c,d,e = (0,0,0,0,0)
    ## matrices
    if 1 == a
        loadPath = script_dir*"\\Data\\CPmats\\M7.txt"
        @testset "loadMatfromtxt and saveMat2txt" begin
            M7 = loadMatfromtxt(loadPath)

            M = M7
            savePath = script_dir*"\\Output\\"
            saveName = "test-example"
            saveMat2txt(M, savePath, saveName )
        end

        @testset "testNN" begin
            M = rand(10,10)
            @test testNN(M) == true
        end

        @testset "testPSD" begin
            a = rand(10,1)
            M = a* tra(a)
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
            A, a_ℓ= genRandCPmatrix(8,5)
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


            MonBase = make_mon_expo(5,4)
            α = [0  0 0 2 2]
            @test  get_mon_index(MonBase, α) == 124
        end

        @testset "make_mom__expo_mat_dict(n::Int,t::Int)" begin

            for n ∈ 6:10 , t ∈ 2:4
                local MonExp = make_mom_expo_mat_dict(n,t)
                local MonBase = make_mon_expo(n,t)
                for key in keys(MonExp)
                    for val in MonExp[key]
                        @test   key == MonBase[val[1]] + MonBase[val[2]]
                    end
                end
            end
        end
    end

    ## This code checks if the constraints are implemented as intended.
    if 1 == c
        A = loadMatfromtxt(loadPath)
        t = 2
        n = size(A)[1]
        Lx = make_dummy_var(n,t)
        che = (n,t) -> binomial(n+1,t)

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

    ## Test Computeξₜᶜᵖ for example matrices.
    if 1 == d
        @testset "Computeξₜᶜᵖ" begin
            cp_mats = ["M11tilde.txt"  "M6.txt"  "M7.txt"  "M7tilde.txt"  "M8tilde.txt"  "M9tilde.txt"]
            loadPath = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\CPmats\\"*cp_mats[3]
            A = loadMatfromtxt(loadPath)

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
end





# A = [1.0 0.5 ; 0.5 1.0]
#
# mom_mat_Temp       = rec_mom_mat(A,t,Lx_1)
#
# Lx_0,model_ξ₂_0ᶜᵖ  = Computeξₜᶜᵖ(A, t, false,0,false)
# objective_value(model_ξ₂_0ᶜᵖ)
# mom_mat_Temp            = rec_mom_mat(A,t,Lx_0)
# eigvals(mom_mat_Temp)
#
#
# Lx_2,model_ξ₂_2ᶜᵖ  = Computeξₜᶜᵖ(A, t, false,2.2,false)
# objective_value(model_ξ₂_2ᶜᵖ)
## The G constriants make_G_con and make_G_con2 are the same upto permutation

# Lx_Temp,model_ξ₂Tempᶜᵖ  = Computeξₜᶜᵖ(A, t, false,0,false)
# mom_mat_Temp            = rec_mom_mat(A,t,Lx_Temp)
# isposdef(mom_mat_Temp)
# println([j for j in eigvals(mom_mat_Temp) if j < 0])



# ξ₂Gᶜᵖ = objective_value(model_ξ₂Gᶜᵖ)
# LxG2,model_ξ₂Gᶜᵖ2  = Computeξₜᶜᵖ(A, t, dag_con,2.2,xx_con)
# ξ₂Gᶜᵖ2 = objective_value(model_ξ₂Gᶜᵖ2)
#
#
#
# LxwG,model_ξ₂wGᶜᵖ  = Computeξₜᶜᵖ(A, t, dag_con,1,xx_con)
# ξ₂wGᶜᵖ = objective_value(model_ξ₂wGᶜᵖ)
# ξ₂wGᶜᵖ = objective_value(model_ξ₂wGᶜᵖ)



# model = Model(Mosek.Optimizer)
# # A = ones(2,2)
# n = size(A)[1]
# t = 2
# list_of_keys = make_mom_expo_keys(n, t)
# @variable(model, Lx[list_of_keys])
#
# G_con                 = make_G_con(A,t,Lx)
# G_con2                = make_G_con2(A,t,Lx)
#
# p = gen_tens_perm(n)
#
# G_con ==  G_con2



##
# value.(LxG) == value.(LxwG)
# mom_matG        = rec_mom_mat(A,t,LxG)
# mom_matwG       = rec_mom_mat(A,t,LxwG)
# mom_matG == mom_matwG
# ##
# mom_matG[2:(1+n),2:(1+n)]  ≈ A
# mom_matwG[2:(1+n),2:(1+n)]  ≈ A
# ## Test the soulution
# # @show eigvals(mom_matG)
# # @show eigvals(mom_matwG)
# #
# @show [i  for i  in eigvals(mom_matwG) if i <  0]
# @show [i  for i  in eigvals(mom_matG) if i <  0]
#
#
# Lxf,model_ξ₂fᶜᵖ  = Computeξₜᶜᵖ(A, t, true,2,true)
# ξ₂fᶜᵖ = objective_value(model_ξ₂fᶜᵖ)
#
#
# function eval_keys(constraint)
#     nonPSD_keys = []
#     for key in keys(loc_con)
#         if isposdef(value.(loc_con[(key)]))
#             append!(nonPSD_keys,[key])
#         end
#     end
#     return nonPSD_keys
# end
#
# loc_con       = make_loc_con(A,t,Lxf)
# nPSD_loc_keys = eval_keys(loc_con)
#
# xx_con        = make_xx_con(A,t,Lxf)
# nPSD_xx_keys  = eval_keys(xx_con)
#
# dag_con       = make_dag_con(A,t,Lxf)
# nPSD_xx_keys  = eval_keys(dag_con)
#
# weakG_con     = make_weakG_con(A,t,Lxf)
# isposdef(weakG_con[1])
#
# G_con         = make_G_con(A,t,Lxf)
# isposdef(G_con)

# @show eigvals(mom_matG) .> eigvals(mom_matwG)
#
# isposdef(mom_matG)
# isposdef(mom_matwG)
# ## Weak Constraints test
# # M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ) ⪰ 0, ℓ ∈ 0,1,t-deg(g)/2
# weakG_conG  =  make_weakG_con(A,t,LxG)
# wG_cG = value.(weakG_conG[1])
# weakG_conwG =  make_weakG_con(A,t,LxwG)
# wG_cwG  = value.(weakG_conwG[1])
#
# G_conwG =  make_G_con(A,t,LxwG)
# G_cwG   = value.(G_conwG)
# G_conG  =  make_G_con(A,t,LxG)
# G_cG    = value.(G_conwG)
#
#
# isposdef(G_cG)
# isposdef(wG_cwG)
# isposdef(wG_cG)
# isposdef(G_cwG)
#
#
# @show eigvals(G_cG)[1:5]
# @show eigvals(wG_cwG)[1:5]
# @show eigvals(wG_cG)[1:5]
# @show eigvals(G_cwG)[1:5]
#
#
#
# G_conwG =  make_G_con(A,t,LxwG)
# G_conG =  make_G_con(A,t,LxG)
# isposdef(value.(G_conG))
#
#
# nar = make_loc_con(A,t,LxG)
#
# for key in keys(nar)
#     if ~isposdef(value.(nar[key]))
#         println(key)
#     end
# end
#
# isposdef(value.(nar[(5, 8)]))
#  eigvals(value.(nar[(5, 8)]))
#
#
# xx_con = make_xx_con(A,t,LxG)
# for key in keys(xx_con)
#     isposdef()
