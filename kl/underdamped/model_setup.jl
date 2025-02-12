using InfiniteOpt, Distributions, Ipopt;
#using Trapz;
using CSV;
using DataFrames;
using ForwardDiff;

#get parameters and boundary conditions
include("../params.jl")

##################################
#Get the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#momentum
@infinite_parameter(model, p in [-10, 10], num_supports = num_supports_p)

#position
@infinite_parameter(model, q in [-3.5, 3.5], num_supports = num_supports_q)


#let's define some good initial guesses
function u_init(t,q)
    return ((q-1)^3)
end
    
function v_init(t,p,q)
     return ((1/4)*((q^2 - 1)^2))
end
    
function rho_init(t,p,q)
    return underdamped_p_initial(p,q)
end


#define the variables
#first remove any existing definitions
#unregister(model, :u)
#unregister(model, :v)
#unregister(model, :rho)


#the optimal control
@variable(model, u, Infinite(t,q), start = u_init)

#the density
@variable(model, rho>=0, Infinite(t,p,q), start = rho_init)

#the optimal control
@variable(model, v, Infinite(t,p,q), start = v_init)

#define the objective
@objective(model, Min, integral(integral(integral((u^2)*rho,t),q), p)/4)

#fokker planck
@constraint(model, deriv(rho,t) + epsilon*p*deriv(rho,q) - deriv((p+u)*rho,q) + deriv(rho,p,p)== 0)

#hjb
@constraint(model, deriv(v,t) - p*deriv(v,p) + deriv(v,p,p) + epsilon*p*deriv(v,q) - epsilon*u*deriv(v,p) == -((epsilon*u)^2)/4)

#stationarity condition
@constraint(model, integral((rho/integral(rho,p))*deriv(v,q),p) == epsilon*u/2)

#boundary conditions on rho
@constraint(model, rho(0,p,q) == underdamped_p_initial(p,q))
@constraint(model, rho(T,p,q) == underdamped_p_final(p,q))

#add constraint to make distribution compact
#@constraint(model, integral(rho,p,q) == 1)
@constraint(model, rho(t,-10,q) == 0)
@constraint(model, rho(t,10,q) == 0)
@constraint(model, rho(t,p,-3.5) == 0)
@constraint(model, rho(t,p,3.5) == 0)


# SOLVE THE MODEL

optimize!(model)

print(termination_status(model))

tgrid = [supports(rho)[i][1] for i in 1:length(supports(rho))]
qgrid = [supports(rho)[i][2] for i in 1:length(supports(rho))];
pgrid = [supports(rho)[i][3] for i in 1:length(supports(rho))];

rhovals = reshape(value(rho),(num_supports_t,num_supports_p,num_supports_q));
vvals = reshape(value(v),(num_supports_t,num_supports_p,num_supports_q));
#uvals = repeat(reshape(value(u),(num_supports_t,1,num_supports_q)),1,num_supports_p,1);

qax = vec(unique(qgrid))
pax = vec(unique(pgrid))
times_vec = vec(unique(tgrid))

#save results to csv
####save a csv file  
file_name = "infiniteopt/ipopt_underdamped_kl_v1.csv"

# Define the header as an array of strings
row = ["t" "p" "q" "v" "rho"]
header = DataFrame(row,["t", "p", "q", "v", "rho"])

# Write the header to a new CSV file
CSV.write(file_name, header;header =false)


for j in 1:num_supports_t
     df = DataFrame([times_vec[j] .*vec(ones(num_supports_p*num_supports_q)),
                    vec(pax),
                    vec(qax),
                    #vec(uvals[j,:,:]),
                    vec(vvals[j,:,:]),
                    vec(rhovals[j,:,:])],
                    ["t", "p", "q", "v", "rho"])

    CSV.write(file_name, df, append =true)
    
end
