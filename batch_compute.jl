#"""save computations to a .txt file """
    using Dates # used in the output file naming
    main_dir = @__DIR__
    include("compute.jl")
    # include("mat_IO.jl")
    include("gen_mats.jl")
    n2d4_map = A -> floor(size(A)[1]^2 /4)

#timestamp = replace(string(now()),":"=>"-")[1:16]
#dataDir = main_dir*"\\Data\\"
#outputDir = check_save_dir(main_dir*"\\", "\\Output\\"*timestamp)




function append_to_md(path,text)
    open( path, "a") do f
        write(f,text)
    end
end

function batchCompξ₂ᶜᵖ(loadName)
    counter = 1
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

        loadPath = loadDir*matfile
        A        = loadMatfromtxt(loadPath)
        n2d4     = n2d4_map(A)
        MatName  = cut_ext(matfile)

        t = 2
        # First three columns
        append_to_md(save_path,"|"*MatName*"|___ |$n2d4|")

        "|Matrix| cp-rank| n²/4|    ξ₂ᶜᵖ|ξ₂wGᶜᵖ| ξ₂Gᶜᵖ|    ξ₂ₓₓᶜᵖ|ξ₂ₓₓwGᶜᵖ|ξ₂ₓₓGᶜᵖ|   ξ₂ᵩᶜᵖ|ξ₂ᵩwGᶜᵖ|ξ₂ᵩGᶜᵖ|   ξ₂ₓₓᵩᶜᵖ|ξ₂ₓₓᵩwGᶜᵖ|ξ₂ₓₓᵩGᶜᵖ| "

        for dag_con in [false,true]
            for xx_con in [false,true]
                for G_con in [0,1,2.1]
                    if dag_con ||  xx_con # || G_con == 0 ||  G_con == 1
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

        println("############################  $counter #############################  ")
        counter = counter + 1
    end

end


function batchModelξ₂ᶜᵖ()
    counter = 1
    for MatName in ["M6", "M7", "M7t", "M8t", "M9t","M11t"]
        A = gen_cp_mats(MatName)
        t = 2
        # First three columns
        for dag_con in [false,true]
            for xx_con in [false,true]
                for G_con in  [0,1,2.1]
                    MatName_T = MatName
                    if xx_con
                        MatName_T = MatName_T*"-xx"
                    end
                    if dag_con
                        MatName_T = MatName_T*"-dag"
                    end
                    if G_con == 1
                        MatName_T = MatName_T*"-wG"
                    end
                    if G_con == 2.1
                        MatName_T = MatName_T*"-G"
                    end
                    model = Modelξₜᶜᵖ(A,t,dag_con,G_con,xx_con)
                    JuMP.write_to_file(model, "DAT-s/$MatName_T.dat-s", format=MOI.FileFormats.FORMAT_SDPA)
                    #println(MatName_T)
                end
            end
        end
        println("############################  $counter #############################  ")
        counter = counter + 1
    end
end

batchModelξ₂ᶜᵖ()


A = gen_cp_mats("M11t")
dataDir_lst =  ["CPmats", "randCPmats", "DDSNNmats","SNNMats"]
for i in 1:0
    batchCompξ₂ᶜᵖ(dataDir_lst[i])
end


# model = Modelξₜᶜᵖ(A,2,false,2.1,false)
#
#
# # Pkg.add("SDPA")
# using SDPA
# set_optimizer(model, SDPA.Optimizer)
# optimize!(model)
#
# println("Primal: ", primal_status(model))
# println("Dual: ", dual_status(model))
# println("Objective: ", objective_value(model))

""" Not used """
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
