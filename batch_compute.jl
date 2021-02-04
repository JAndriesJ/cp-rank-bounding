#"""save computations to a .txt file """
using Dates # used in the output file naming
main_dir = @__DIR__
include("compute.jl")
include("mat_IO.jl")

timestamp = replace(string(now()),":"=>"-")[1:16]

dataDir = main_dir*"\\Data\\"
outputDir = check_save_dir(main_dir*"\\", "\\Output\\"*timestamp)
n2d4_map = A -> floor(size(A)[1]^2 /4)



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

        println("############################  $counter #############################  ")
        counter = counter + 1
    end

end


dataDir_lst =  ["CPmats", "randCPmats", "DDSNNmats","SNNMats"]
for i in 1:1
    batchCompξ₂ᶜᵖ(dataDir_lst[i])
end







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
