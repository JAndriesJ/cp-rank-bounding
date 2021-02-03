#"""save computations to a .txt file """
using Dates # used in the output file naming
    main_dir = @__DIR__
    include("compute.jl")
    include("mat_IO.jl")

    timestamp = replace(string(now()),":"=>"-")[1:16]

    dataDir = main_dir*"\\Data\\"
    outputDir = check_save_dir(main_dir*"\\", "\\Output\\"*timestamp)


function append_to_md(path,text)
    open( path, "a") do f
        write(f,text)
    end
end

n2d4_map = A -> floor(size(A)[1]^2 /4)

function batchCompξ₂ᶜᵖ(loadName)
    counter = 1
    #colnames  = "Matrix,    cp-rank,    n²/4,  ξ₂ᶜᵖ,      ξ₂ᵩᶜᵖ,      ξ₂ᵩweak⊗ᶜᵖ,     ξ₂ᵩ⊗ᶜᵖ,    ξ₂ᵩ⊗ᶜᵖ + xᵢxⱼ"
    Header1 = "|Matrix| cp-rank| n²/4|    ξ₂ᶜᵖ|ξ₂wGᶜᵖ| ξ₂Gᶜᵖ|    ξ₂ₓₓᶜᵖ|ξ₂ₓₓwGᶜᵖ|ξ₂ₓₓGᶜᵖ|   ξ₂ᵩᶜᵖ|ξ₂ᵩwGᶜᵖ|ξ₂ᵩGᶜᵖ|   ξ₂ₓₓᵩᶜᵖ|ξ₂ₓₓᵩwGᶜᵖ|ξ₂ₓₓᵩGᶜᵖ| "
    Header2 = "|---   |---     |---  |---     |---   |---   |---       |---     |---    |---     |---    |---   |---       |---      |---     |"
    loadDir = dataDir*loadName*"\\"
    Matfiles =  cd(readdir, loadDir)
    open( outputDir*loadName*".md", "w") do f
        write(f, Header1 *  "  \n")
        write(f, Header2 *  "  \n")
    end
    save_path = outputDir*loadName*".md"
    for matfile in Matfiles
        # if counter > 3
        #    break
        # end
        loadPath = loadDir*matfile
        A        = loadMatfromtxt(loadPath)
        n2d4     = n2d4_map(A)
        MatName  = cut_ext(matfile)


        #ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂weakTensᶜᵖ,ξ₂ᵩweakTensᶜᵖ,ξ₂Tensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
        #ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂weakTensᶜᵖ,ξ₂Tensᶜᵖ = compξᶜᵖ_1(M,t)
        t = 2
        # First three columns
        append_to_md(save_path,"|"*MatName*"|___ |$n2d4|")

        "|Matrix| cp-rank| n²/4|    ξ₂ᶜᵖ|ξ₂wGᶜᵖ| ξ₂Gᶜᵖ|    ξ₂ₓₓᶜᵖ|ξ₂ₓₓwGᶜᵖ|ξ₂ₓₓGᶜᵖ|   ξ₂ᵩᶜᵖ|ξ₂ᵩwGᶜᵖ|ξ₂ᵩGᶜᵖ|   ξ₂ₓₓᵩᶜᵖ|ξ₂ₓₓᵩwGᶜᵖ|ξ₂ₓₓᵩGᶜᵖ| "

        for dag_con in [false,true]
            for xx_con in [false,true]
                for G_con in [0,1,2.2]
                    if dag_con ||  xx_con  || G_con == 0 ||  G_con == 1
                        append_to_md(save_path,"skip|")
                     else
                        Lx,model_ξ₂ᶜᵖ = Computeξₜᶜᵖ(A, t, dag_con, G_con, xx_con)
                        ξₜᶜᵖvar       = round(objective_value(model_ξ₂ᶜᵖ), digits = 4)
                        append_to_md(save_path,"$ξₜᶜᵖvar|")
                    end
                end
            end
        end
        append_to_md(save_path," \n")



        # if ~(ξ₂ᶜᵖ <= ξ₂ᵩᶜᵖ<= ξ₂ᵩweakTensᶜᵖ<= ξ₂ᵩTensᶜᵖ<= ξ₂ᵩTensₓₓᶜᵖ)
        #     ord_vio = "Order-Violation-"
        # end

        # open( outputDir*loadName*".md", "a") do f
        #     write(f,line)
        # end


        println("############################  $counter #############################  ")
        counter = counter + 1
    end

end


dataDir_lst =  ["CPmats", "randCPmats", "DDSNNmats","SNNMats"]
for i in 1:1
    batchCompξ₂ᶜᵖ(dataDir_lst[i])
end

# """TO DO: make a batch update """
# replace(string(now()),":"=>"-")


function batchCompStuff(loadName)
    counter = 1
    loadDir = dataDir*loadName*"\\"
    Matfiles =  cd(readdir, loadDir)
    for matfile in Matfiles
        if counter > 4
           break
        end
        loadPath = loadDir*matfile
        M       = loadMatfromtxt(loadPath)

        if  ~testNN(M) || ~testPSD(M)
            println(matfile)
        end


        counter = counter + 1
    end
end



#
# # 4th column
# Lx,model_ξ₂ᶜᵖ           = Computeξₜᶜᵖ(A, t, false,0,false)
# append_to_md(save_path,"$ξ₂ᶜᵖ|")
# # 5th column
# Lx,model_ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(M, t, true, 0,false)
# append_to_md(save_path,"$ξ₂ᵩᶜᵖ|")
# # 6th column
# Lx,model_ξ₂ₓₓᶜᵖ         = Computeξₜᶜᵖ(M, t, false, 0, true)
# append_to_md(save_path,"$ξ₂ₓₓᶜᵖ|")
# # 7th column
# Lx,model_ξ₂weakGᶜᵖ   = Computeξₜᶜᵖ(M, t, false, 1,false)
# append_to_md(save_path,"$ξ₂weakGᶜᵖ|")
# # 8th column
# Lx,model_ξ₂ᵩweakGᶜᵖ  = Computeξₜᶜᵖ(M, t, true, 1,false)
# append_to_md(save_path,"$ξ₂ᵩweakGᶜᵖ|")
# # 9th column
# Lx,model_ξ₂Gᶜᵖ       = Computeξₜᶜᵖ(M, t, false, 2,false)
# append_to_md(save_path,"$ξ₂Gᶜᵖ|")
# # 10th column
# Lx,model_ξ₂ᵩGᶜᵖ      = Computeξₜᶜᵖ(M, t, true, 2,false)
# append_to_md(save_path,"$ξ₂ᵩGᶜᵖ|")
# # 11th column
# Lx,model_ξ₂ᵩGₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, true, 2, true)
# ξ₂ᵩGₓₓᶜᵖ
# append_to_md(save_path,"$ξ₂ᵩGₓₓᶜᵖ| \n")
