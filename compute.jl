using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.
# using GLPK
# using Cbc
include("constraints.jl")

"""This is where the ξₜᶜᵖ is calculated for matrix A """
function Computeξₜᶜᵖ(A,t,isDag,GtensL,isXX)
    n = size(A)[1]
    ## Begin making the model
    model = Model(Mosek.Optimizer)
    # model = Model(GLPK.Optimizer)
    # Define all variables that occur in the moment matrix.
    list_of_keys = make_mom_expo_keys(n, t)
    @variable(model, Lx[list_of_keys] )
## Build the moment matrix and constrain it to be PSD.
    mom_matₜ_expo = make_mon_expo_mat(n,t,true)
    mom_matₜ      = index_to_var(Lx, mom_matₜ_expo)
    @SDconstraint(model, mom_con, mom_matₜ >= zeros(size(mom_matₜ)))
## Second order Moment constraints
    mom_mat₌₁ = make_mon_expo_mat(n,1,false)
    Lx_mom_mat₌₁ = index_to_var(Lx,mom_mat₌₁)
    @constraint(model, fix_con,Lx_mom_mat₌₁ .==  A)

    # for i ∈ 1:n, j ∈ 1:n
    #     eᵢ = get_std_base_vec(n,i)
    #     eⱼ = get_std_base_vec(n,j)
    #     @constraint(model, Lx[eᵢ + eⱼ] ==  A[i,j])
    #     #@constraint(model, Lx_mom_mat₌₁[i,j] ==  A[i,j])
    # end
## Localizing g constraint
    loc_con = make_loc_con(A,t,Lx)
    Z_mat = zeros(size(loc_con[(1,1)]))
    for key in keys(loc_con)
        @SDconstraint(model, loc_con[key] >= Z_mat)
    end
## Dagger constraints
    if isDag
        println("----------------Dagger constraints are active")
        dag_con  = make_dag_con(A,t,Lx)
        for key in keys(dag_con)
            @constraint(model, dag_con[key] .>= zeros(size(dag_con[key])))
        end
    end

## Localizing XX constraint
    if isXX
        println("----------------XX constraints are active")
        xx_con = make_xx_con(A,t,Lx)
        for k in 1:n
            for h in (k+1):n
                @SDconstraint(model, xx_con[(k,h)] >= zeros(size(xx_con[(k,h)])))
            end
        end
    end
## G Constraints

    if GtensL == 1
        println("----------------Weak G-constraints are active")
        weakG_con = make_weakG_con(A,t,Lx)
        for key in keys(weakG_con)
            weakG_con_key = weakG_con[key]
           @SDconstraint(model, weakG_con_key >= zeros(size(weakG_con_key)))
        end
    end
    if GtensL == 2
        println("----------------G-constraints are active")
        G_con                 = make_G_con(A,t,Lx)
        @SDconstraint(model, G_con >= zeros(size(G_con)))
    end

    #  Set objective
    @objective(model, Min, Lx[zeros(n)] )
    # optimize
    optimize!(model)

    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    # println("Variables: ", "x = ",value(x))
    if string(primal_status(model)) == "NO_SOLUTION"
        return NaN
    else
        # return
        return Lx,model
    end
     # termination_status(model)
     # value.(Lx)
     # dual.(model)
     # round(objective_value(model), digits = 8)
end


function rec_mom_mat(n::Int64,t::Int64,Lx)
    MB_exp  = make_mon_expo_mat(n,t)
    MB      = index_to_var(Lx, MB_exp)
    mom_mat = value.(MB)
    return mom_mat
end

function rec_mom_mat(A::Array{Float64,2},t::Int64,Lx)
    @assert size(A)[1] == size(A)[2]
    n = size(A)[1]
    return rec_mom_mat(n::Int64,t::Int64,Lx)
end
