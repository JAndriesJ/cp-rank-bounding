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


function batchCompξ₂ᶜᵖ(loadName)
    counter = 1
    #colnames  = "Matrix,    cp-rank,    n²/4,  ξ₂ᶜᵖ,      ξ₂ᵩᶜᵖ,      ξ₂ᵩweak⊗ᶜᵖ,     ξ₂ᵩ⊗ᶜᵖ,    ξ₂ᵩ⊗ᶜᵖ + xᵢxⱼ"
    Header1 = "|Matrix| cp-rank| n²/4| ξ₂ᶜᵖ|  ξ₂ᵩᶜᵖ| ξ₂ₓₓᶜᵖ | ξ₂weakGᶜᵖ|ξ₂ᵩweakGᶜᵖ| ξ₂Gᶜᵖ|ξ₂ᵩGᶜᵖ| ξ₂ᵩGᶜᵖ + xᵢxⱼ|"
    Header2 = "|---|---|---|---|---|---|---|---|---|---|---|"
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
        M        = loadMatfromtxt(loadPath)
        n2d4     = A -> floor(size(A)[1]^2 /4)
        MatName  = cut_ext(matfile)


        #ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂weakTensᶜᵖ,ξ₂ᵩweakTensᶜᵖ,ξ₂Tensᶜᵖ,ξ₂ᵩTensᶜᵖ, ξ₂ᵩTensₓₓᶜᵖ = compξᶜᵖ(M,2)
        #ξ₂ᶜᵖ, ξ₂ᵩᶜᵖ, ξ₂weakTensᶜᵖ,ξ₂Tensᶜᵖ = compξᶜᵖ_1(M,t)
        t = 2

        append_to_md(save_path,"|"*MatName*"|___ |$n2d4|")
        ξ₂ᶜᵖ           = Computeξₜᶜᵖ(M, t, false,0,false)
        append_to_md(save_path,"$ξ₂ᶜᵖ|")
        ξ₂ᵩᶜᵖ          = Computeξₜᶜᵖ(M, t, true, 0,false)
        append_to_md(save_path,"$ξ₂ᵩᶜᵖ|")
        ξ₂ₓₓᶜᵖ         = Computeξₜᶜᵖ(M, t, false, 0, true)
        append_to_md(save_path,"$ξ₂ₓₓᶜᵖ|")
        ξ₂weakGᶜᵖ   = Computeξₜᶜᵖ(M, t, false, 1,false)
        append_to_md(save_path,"$ξ₂weakGᶜᵖ|")
        ξ₂ᵩweakGᶜᵖ  = Computeξₜᶜᵖ(M, t, true, 1,false)
        append_to_md(save_path,"$ξ₂ᵩweakGᶜᵖ|")
        ξ₂Gᶜᵖ       = Computeξₜᶜᵖ(M, t, false, 2,false)
        append_to_md(save_path,"$ξ₂Gᶜᵖ|")
        ξ₂ᵩGᶜᵖ      = Computeξₜᶜᵖ(M, t, true, 2,false)
        append_to_md(save_path,"$ξ₂ᵩGᶜᵖ|")
        ξ₂ᵩGₓₓᶜᵖ    = Computeξₜᶜᵖ(M, t, true, 2, true)
        append_to_md(save_path,"$ξ₂ᵩGₓₓᶜᵖ| \n")


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

# loadpath = "C:\\Users\\andries\\all-my-codes\\cp-rank-bounding\\Output\\CPmats.md"
# open( loadpath, "a") do f
#     write(f,"11111111111111 \n")
# end
#






dataDir_lst =  ["CPmats", "randCPmats", "DDSNNmats","SNNMats"]
for i in 1:4
    batchCompξ₂ᶜᵖ(dataDir_lst[i])
end

# """TO DO: make a batch update """
# replace(string(now()),":"=>"-")
