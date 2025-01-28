using InfiniteOpt, Distributions, Ipopt;
using Trapz;
using CSV;
using DataFrames;
using ForwardDiff;


#calculate normalisation for distributions
norm_range = Array(range(-8,8,8000))

global norminitial = abs.(trapz(norm_range,exp.(-((norm_range.-1).^4)/4)))*sqrt(2*pi)
global normfinal = abs.(trapz(norm_range,exp.(-(((norm_range.^2).-1).^2)/4)))*sqrt(2*pi)

#####boundary conditions#######
global function p_initial(p,q)

    return exp.(-((q.-1).^4)/4) .* exp.(-(p.^2)/2) /norminitial
end;


global function p_final(p,q)

    return exp.(-(((q.^2).-1).^2)/4) .* exp.(-(p.^2)/2) /normfinal
end;

##LATTICE PARAMETERS!!!
T = 0.2
num_supports_t = 11
num_supports_q = 10
num_supports_p = 10

##################################
#Get the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#momentum
@infinite_parameter(model, p in [-10, 10], num_supports = num_supports_p)

#position
@infinite_parameter(model, q in [-3, 3], num_supports = num_supports_q)


#let's define some good initial guesses
function u_init(t,q)
    return ((q-1)^3)
end
    
function v_init(t,q)
     return ((1/4)*((y^2 - 1)^2))
end
    
function rho_init(t,p,q)
    return p_initial(p,q)
end


#define the variables
#first remove any existing definitions
#unregister(model, :u)
#unregister(model, :v)
#unregister(model, :rho)


#the optimal control
@variable(model, u, Infinite(t,p,q), start = u_init)

#the density
@variable(model, rho>=0, Infinite(t,p,q), start = rho_init)

#the optimal control
@variable(model, v, Infinite(t,p,q), start = v_init)

#define the objective
@objective(model, Min, integral(integral((u^2)*rho,q), t)/4)

#fokker planck
@constraint(model, deriv(rho,t) - deriv(u*rho,q) - deriv(rho,q,q)== 0)

#hjb
@constraint(model, deriv(v,t) - u*deriv(v,q) + deriv(v,q,q) + (u^2)/4 == 0)

#stationarity condition
@constraint(model, integral((rho(t,p,q)/integral(rho(t,p,q),p))*deriv(v,q),p) == u)

#boundary conditions on rho
@constraint(model, rho(0,p,q) == p_initial(p,q))
@constraint(model, rho(T,p,q) == p_final(p,q))

#add constraint to make distribution compact
#@constraint(model, integral(rho,p,q) == 1)
@constraint(model, rho(t,-10,q) == 0)
@constraint(model, rho(t,10,q) == 0)
@constraint(model, rho(t,p,-1) == 0)
@constraint(model, rho(t,p,3) == 0)


# SOLVE THE MODEL
optimize!(model)

print(termination_status(model))

tgrid = [supports(rho)[i][1] for i in 1:length(supports(rho))]
qgrid = [supports(rho)[i][2] for i in 1:length(supports(rho))];
pgrid = [supports(rho)[i][3] for i in 1:length(supports(rho))];

rhovals = reshape(value(rho),(num_supports_t,num_supports_p,num_supports_q));
vvals = reshape(value(v),(num_supports_t,num_supports_p,num_supports_q));
uvals = reshape(value(u),(num_supports_t,num_supports_p,num_supports_q));

qax = vec(unique(qgrid))
pax = vec(unique(pgrid))
times_vec = vec(unique(tgrid))

#save results to csv
####save a csv file  
file_name = "infiniteopt/ipopt_underdampedkl_v1.csv"

# Define the header as an array of strings
row = ["t" "x" "du" "v" "rho"]
header = DataFrame(row,["t", "x", "du", "v", "rho"])

# Write the header to a new CSV file
CSV.write(file_name, header;header =false)


for j in 1:num_supports_t
     df = DataFrame([times_vec[j] .*vec(ones(num_supports_p*num_supports_q)),
                    vec(qax),
                    vec(uvals[j,:,:]),
                    vec(vvals[j,:,:]),
                    vec(rhovals[j,:,:])],
                    ["t", "x", "du", "v", "rho"])

    CSV.write(file_name, df, append =true)
    
end
