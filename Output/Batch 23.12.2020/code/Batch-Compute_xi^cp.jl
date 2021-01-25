#"""save computations to a .txt file """
include("Compute_xi^cp.jl")
include("MatUtils.jl")

dataDir = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\"
outputDir = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Output\\"





function batchCompξ₂ᶜᵖ(loadName,SaveName)
    counter = 1
    colnames  = "Matrix,    cp-rank,    n²/4,  ξ₂ᶜᵖ,      ξ₂ᵩᶜᵖ,      ξ₂ᵩweak⊗ᶜᵖ,     ξ₂ᵩ⊗ᶜᵖ,    ξ₂ᵩ⊗ᶜᵖ + xᵢxⱼ"
    loadDir = dataDir*loadName
    Matfiles =  cd(readdir, loadDir)
    file_name = SaveName
    open( outputDir*file_name*".txt", "w") do f
        write(f, colnames *  "  \n")
        for matfile in Matfiles
            loadPath = loadDir*matfile
            M = loadMatfromtxt(loadPath)
            n2d4          = n2div4(M)
            MatName       = extName(matfile)
            println(MatName )
            ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)

            if any(ξ₂ᶜᵖ <= ξ₂ᵩᶜᵖ<= ξ₂ᵩweakTensᶜᵖ<= ξ₂ᵩTensᶜᵖ<= ξ₂ᵩTensₓₓᶜᵖ)
                MatName = "error-"*MatName
            end
            write(f, MatName*", ___     ,$n2d4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
            println("############################  $counter #############################  ")
            counter = counter + 1
        end
    end
end

# batchCompξ₂ᶜᵖ("CPmats\\","CPmats")
batchCompξ₂ᶜᵖ("CPrandMat"*"\\","CPrandMats")
batchCompξ₂ᶜᵖ("DDSNN_Mat"*"\\","DDSNN_Mat")
batchCompξ₂ᶜᵖ("SNN_Mat"*"\\","SNN_Mat")
#








# if false
#     #cp_ranks = [6 14 14 18 26 34]
#     # cp_matrices = [M6, M7, M7tilde, M8tilde, M9tilde, M11tilde]
#     # cp_matrices_names = ["M6", "M7", "M7tilde", "M8tilde", "M9tilde", "M11tilde"]
#
#     loadDir = dataDir*"CPmats\\"
#     Matfiles =  cd(readdir, loadDir)
#     file_name = "CPmats"
#     open( outputDir*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for matfile in Matfiles
#             loadPath = loadDir*matfile
#             M = loadMatfromtxt(loadPath)
#
#             n²div4        = n²div4(cp_matrices[ind])
#             cp_ranks_text = string(cp_ranks[ind])
#             MatName       = extName(matfile)
#             ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(cp_matrices[ind],2)
#
#             write(f,  MatName*","* cp_ranks_text*"     ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#             println("############################  $counter #############################  ")
#             counter = counter + 1
#         end
#     end
# end



# # Explotration 1: MORE random matrices
# if false
#     open( SavePath *"_Random_NN_outer_prod_"*file_name*".txt", "w") do f
#         # Explotration 1: MORE random matrices
#         write(f, colnames *  "  \n")
#         for  n in [6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9]
#             r = rand(4:11)
#             M, a_ℓ = genCPmatrix(n,r)
#             n²div4 =  floor(n^2 / 4)
#             Mtext  = "$n Rᶜᵖ_$r"
#             ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#             write(f,  Mtext*" , $r, ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#             println("############################  $counter #############################  ")
#             counter = counter + 1
#         end
#     end
# end
# # Explotration 2: For real, where the fat matrices at?
# if false
#     open( SavePath *"Random_NN_outer_prod_high_rank"*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for  n in [6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9]
#             r = rand(30:40)
#             M, a_ℓ = genCPmatrix(n,r)
#             n²div4 =  floor(n^2 / 4)
#             Mtext  = "$n Rᶜᵖ_$r"
#             ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#             write(f,  Mtext*" , $r, ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#             println("############################  $counter #############################  ")
#             counter = counter + 1
#         end
#     end
# end
#
# # Explotration 3: MORE random matrices
# if false
#     open( SavePath *"_Random_tensor_prod_construction"*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for  n in [2,2,2,2,2,2,2,2,2,2]
#            r = rand(1:5)
#            M, a_ℓ = genVariantCPmatrix(n,r)
#            n²div4 =  floor(n^2 / 4)
#            Mtext = "$n ³VRᶜᵖ_$r"
#            ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#            write(f,  Mtext*" , ?$r?, ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#            println("############################  $counter #############################  ")
#             counter = counter + 1
#        end
#     end
# end
#
# # Explotration 4: What do symetric nonnegatvie with
# if false
#     open( SavePath *"_sym_NN_"*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for  n in [7,7,7,8,8,8,9,9,9]
#            M = genSNNmatrix(n)
#            n²div4 =  floor(n^2 / 4)
#            Mtext = "$n _sym_NN_"
#            ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#            write(f,  Mtext*" , ?? ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#            println("############################  $counter #############################  ")
#             counter = counter + 1
#        end
#     end
# end
#
# # Explotration 4: What do symetric nonnegatvie with
# if false
#     open( SavePath *"_SNN_diag_shift_"*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for  n in [7,8,9]
#            M = genSNNmatrix(n)
#            for λ in 0:0.25:1
#                M = DominateDiagonal(M,λ)
#                n²div4 =  floor(n^2 / 4)
#                Mtext = "$n _sym_NN_λ$λ"
#                ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#                write(f,  Mtext*" , ?? ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#                println("############################  $counter #############################  ")
#                counter = counter + 1
#             end
#        end
#     end
# end
#
# if true
#     open( SavePath *"_Bi_part_graph_"*file_name*".txt", "w") do f
#         write(f, colnames *  "  \n")
#         for  p in [3,3,4,4,4,4,5,5,5]
#            q = rand(3:5)
#            a = rand(1:10)
#            b = rand(1:10)
#            M = BipartiteGraphclass(p,q,a,b)
#            n²div4 =  floor((p+q)^2 / 4)
#            cprankbound = p*q + p + q
#            Mtext = "$p _$q _Bi_part_g_$a _$b"
#            ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂ᵩweakTensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
#            write(f,  Mtext*" , $cprankbound ,$n²div4,       $ξ₂ᶜᵖ,      $ξ₂ᵩᶜᵖ,      $ξ₂ᵩweakTensᶜᵖ,     $ξ₂ᵩTensᶜᵖ,      $ξ₂ᵩTensₓₓᶜᵖ \n")
#
#            println("############################  $counter #############################  ")
#             counter = counter + 1
#        end
#     end
# end




# BatchCompξᶜᵖ("C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Output\\")
