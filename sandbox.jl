using JuMP
using MosekTools

# Needs MatLab
#] add https://github.com/blegat/CDCS.jl.git
#using CDCS

# Need to figure out how to install
# Has special structure
#using CDCS

#Pkg.add("COSMO")
#using COSMO
#Pkg.rm("COSMO")

#Pkg.add("Hypatia")
#using Hypatia
#Pkg.rm("Hypatia")

# Gives StackOverflowError when querying the solution
#Pkg.add("ProxSDP")
#using ProxSDP
#Pkg.rm("ProxSDP")

#Pkg.add("SCS")
#using SCS
#Pkg.rm("SCS")

# Got errors and haven't bothered to get it working
#Pkg.add("SDPAFamily")
#Pkg.rm("SDPAFamily")
#using SDPAFamily


model = Model(ProxSDP.Optimizer)
@variable(model, x)
poes = [2 x; x 1]
@SDconstraint(model, poes >= zeros(2, 2))
@objective(model, Max, x)
optimize!(model)


println("Primal: ", primal_status(model))
println("Dual: ", dual_status(model))
println("Objective: ", objective_value(model))



prog = "1 + 1"
ex1 = Meta.parse(prog)
typeof(ex1)
ex1.head
dump(ex1)
eval(ex1)


prog2 = "Model(Hypatia.Optimizer)"
ex2 = Meta.parse(prog2)
model = eval(ex2)
