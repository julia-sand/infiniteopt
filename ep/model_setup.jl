#get all the packies
using InfiniteOpt, Distributions, Ipopt;
#using Trapz;
using CSV;
using DataFrames;
using ForwardDiff;

include("params.jl")

##################################
#Get the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#space
@infinite_parameter(model, q in [-3.5, 3.5], num_supports = num_supports_q)



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
#@variable(model, u, Infinite(t,q), start = u_init)

#define the objective
#@objective(model, Min, integral(integral(((u^2)-deriv(u,q))*rho,q), t))
@objective(model, Min, integral(integral(rho*(deriv(sigma,q)^2), q),t))

#equation for density 
@constraint(model, deriv(rho,t) - (epsilon^2)*deriv(sigma,q)*deriv(rho,q)== 0)

#burgers equation
@constraint(model, deriv(sigma,t) - ((epsilon^2)/2)*(deriv(sigma,q)^2)== 0)

#stationarity condition
#@constraint(model, deriv(sigma,q) - (deriv(rho,q)/rho) == u)

#boundary conditions on rho
@constraint(model, rho(0,q) == exp(-((q-1)^4)/4)/norminitial)
@constraint(model, rho(T,q) == exp(-((q^2-1)^2)/4)/normfinal)

#normalisation condition
@constraint(model, integral(rho,q) == 1)
@constraint(model, rho(t,-3.5) == 0)
@constraint(model, rho(t,3.5) == 0)


# SOLVE THE MODEL
optimize!(model)


print(termination_status(model))

tgrid = [supports(rho)[i][1] for i in 1:length(supports(rho))]
qgrid = [supports(rho)[i][2] for i in 1:length(supports(rho))];

rhovals = reshape(value(rho),(num_supports_t,num_supports_q));
sigvals = reshape(value(sigma),(num_supports_t,num_supports_q));
#uvals = reshape(value(u),(num_supports_t,num_supports_q));

qax = vec(unique(qgrid))
times_vec = vec(unique(tgrid))

#save results to csv
####save a csv file  
file_name = "infiniteopt/ep/ipopt_overdampedep_v1.csv"

# Define the header as an array of strings
row = ["t" "x" "sigma" "rho"]
header = DataFrame(row,["t", "x", "sigma", "rho"])

# Write the header to a new CSV file
CSV.write(file_name, header;header =false)


for j in 1:num_supports_t
     df = DataFrame([times_vec[j] .*vec(ones(num_supports_q)),
                    vec(qax),
                    #vec(uvals[j,:]),
                    vec(sigvals[j,:]),
                    vec(rhovals[j,:])],
                    ["t", "x", "sigma", "rho"])

    CSV.write(file_name, df, append =true)
    
end



