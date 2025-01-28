#get all the packies
using InfiniteOpt, Distributions, Ipopt;
using Trapz;
using CSV;
using DataFrames;
using ForwardDiff;


####boundary conditions
#intial and final distributions
###########
global function p_initial(y)

    return exp.(-((y.-1).^4)/4)
end;


global function p_final(y)

    return exp.(-(((y.^2).-1).^2)/4);
end;

#calculate normalisation for distributions
norm_range = Array(range(-8,8,8000))

global normfinal = abs.(trapz(norm_range,p_final(norm_range)));
global norminitial = abs.(trapz(norm_range,p_initial(norm_range)));

##LATTICE PARAMETERS!!!
T0 = 2 #0.2
num_supports_t = 11
num_supports_q = 10

epsilon = 0.2
g = 0.01
T = (epsilon^2)*T0

##################################
#Get the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#space
@infinite_parameter(model, q in [-3, 3], num_supports = num_supports_q)



#let's define some good initial guesses
function u_init(t,y)
    return ((y-1)^3)
end
    
function v_init(t,y)

    return ((1/4)*((y^2 - 1)^2)) #/10
end
    
function rho_init(t,y)
    return p_initial(y)/norminitial
end


#define the variables
#first remove any existing definitions
#unregister(model, :u)
#unregister(model, :v)
#unregister(model, :rho)


#the burgers velocity
@variable(model, sigma, Infinite(t,q), start = u_init)

#the density
@variable(model, rho>=0, Infinite(t,q), start = rho_init)

#the optimal control
@variable(model, u, Infinite(t,q), start = u_init)

#define the objective
@objective(model, Min, integral(integral(((u^2)-deriv(u,q))*rho,q), t))

#equation for density 
@constraint(model, deriv(rho,t) - (epsilon^2)*deriv(sigma,q)*deriv(rho,q)== 0)

#burgers equation
@constraint(model, deriv(sigma,t) - ((epsilon^2)/2)*(deriv(sigma,q)^2)== 0)

#stationarity condition
@constraint(model, deriv(sigma,q) - (deriv(rho,q)/rho) == u)

#boundary conditions on rho
@constraint(model, rho(0,q) == exp(-((q-1)^4)/4)/norminitial)
@constraint(model, rho(T,q) == exp(-((q^2-1)^2)/4)/normfinal)

#normalisation condition
@constraint(model, integral(rho,q) == 1)
@constraint(model, rho(t,-3) == 0)
@constraint(model, rho(t,3) == 0)


# SOLVE THE MODEL
optimize!(model)


print(termination_status(model))

tgrid = [supports(rho)[i][1] for i in 1:length(supports(rho))]
qgrid = [supports(rho)[i][2] for i in 1:length(supports(rho))];

rhovals = reshape(value(rho),(num_supports_t,num_supports_q));
sigvals = reshape(value(sigma),(num_supports_t,num_supports_q));
uvals = reshape(value(u),(num_supports_t,num_supports_q));

qax = vec(unique(qgrid))
times_vec = vec(unique(tgrid))

#save results to csv
####save a csv file  
file_name = "infiniteopt/ep/ipopt_overdampedep_v1.csv"

# Define the header as an array of strings
row = ["t" "x" "du" "sigma" "rho"]
header = DataFrame(row,["t", "x", "du", "sigma", "rho"])

# Write the header to a new CSV file
CSV.write(file_name, header;header =false)


for j in 1:num_supports_t
     df = DataFrame([times_vec[j] .*vec(ones(num_supports_q)),
                    vec(qax),
                    vec(uvals[j,:]),
                    vec(sigvals[j,:]),
                    vec(rhovals[j,:])],
                    ["t", "x", "du", "sigma", "rho"])

    CSV.write(file_name, df, append =true)
    
end



