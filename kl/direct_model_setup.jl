#get all the packies
using InfiniteOpt, Distributions, Ipopt;
#using Trapz;
using CSV;
using DataFrames;
using ForwardDiff;

include("params.jl")
#=
Setup the infinite opt model 
=#



#let's define some good initial guesses
function u_init(t,y)
    return ((y-1)^3)
end
    
#function v_init(t,y)
#     return ((1/4)*((y^2 - 1)^2))
#end
    
function rho_init(t,y)
    return p_initial(y)/norminitial
end


##################################
#Setup the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

set_optimizer_attribute(model, "max_iter", 1000)

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#noise
@infinite_parameter(model, xi ~ Normal(0, 1), num_supports = num_supports_q)


@variable(model, q, Infinite(t,xi))

#the optimal control
@variable(model, u, Infinite(t,xi), start = u_init)

#the density
@variable(model, rho>=0, Infinite(t,q), start = rho_init)

#define the objective
@objective(model, Min, integral(integral((u^2)*rho,q), t))

#boundary conditions on rho
@constraint(model, rho(0,q) == p_initial(q)/norminitial)
@constraint(model, rho(T,q) == p_final(q)/normfinal)

#normalisation condition
@constraint(model, integral(rho,q) == 1)
@constraint(model, rho(t,-5) == 0)
@constraint(model, rho(t,5) == 0)


# SOLVE THE MODEL
optimize!(model)


print(termination_status(model))

tgrid = [supports(rho)[i][1] for i in 1:length(supports(rho))]
qgrid = [supports(rho)[i][2] for i in 1:length(supports(rho))];

rhovals = reshape(value(rho),(num_supports_t,num_supports_q));
vvals = reshape(value(v),(num_supports_t,num_supports_q));
uvals = reshape(value(u),(num_supports_t,num_supports_q));

qax = vec(unique(qgrid))
times_vec = vec(unique(tgrid))

#save results to csv
####save a csv file  
file_name = "infiniteopt/ipopt_overdampedkl_v1.csv"

# Define the header as an array of strings
row = ["t" "x" "du" "v" "rho"]
header = DataFrame(row,["t", "x", "du", "v", "rho"])

# Write the header to a new CSV file
CSV.write(file_name, header;header =false)


for j in 1:num_supports_t
     df = DataFrame([times_vec[j] .*vec(ones(num_supports_q)),
                    vec(qax),
                    vec(uvals[j,:]),
                    vec(vvals[j,:]),
                    vec(rhovals[j,:])],
                    ["t", "x", "du", "v", "rho"])

    CSV.write(file_name, df, append =true)
    
end



