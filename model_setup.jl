#get all the packies
using InfiniteOpt, Distributions, Ipopt;
using Plots;
using Trapz;
using CSV;
using DataFrames;
using Interpolations;
using ForwardDiff;
using ImageFiltering;
#using Polynomials;
##filtering functions

#get the normalisation constants




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


#intial and final distributions
###########
global function underdamped_p_initial(p,q)

    return exp.(-(((q.-1).^4)/4 .+ (p.^2)/2))/norminitial
end;


global function underdamped_p_final(p,q)

    return exp.(-((((q.^2).-1).^2)/4 .+ (p.^2)/2))/normfinal;
end;

#get Caluya-Halder output
df_new = CSV.read("results_kl_land_new1.csv",DataFrame,header=true);

time_interval = vec(unique(df_new.t))
global time_interval = round.(time_interval;digits = 4)

global time_grid = Int(length(time_interval));


global time_steps_vec = [round.(time_interval[Int(k+1)]-time_interval[Int(k)];digits =3) for k in 1:length(time_interval)-1];
global plot_times = time_interval[begin:4:end]; #these times are used in the plots

##LATTICE PARAMETERS!!!
T = 1 #0.2
num_supports_t = 11
num_supports_q = 10
num_supports_p = 10

##################################
#Get the model
#create a new model object
model = InfiniteModel(Ipopt.Optimizer);

#time
@infinite_parameter(model, t in [0, T], num_supports = num_supports_t)

#space
@infinite_parameter(model, q in [-3, 3], num_supports = num_supports_q)

#space
@infinite_parameter(model, p in [-3, 3], num_supports = num_supports_p)


#let's define some good initial guesses
function u_init(t,y)
    return ((y-1)^3)/10
end
    
function v_init(t,y)
    #p = Polynomial([0.012501817,1.7759079, 
    #        -0.98531216, 0.29115227, 
    #        0.16276835, -0.04258577,0.01356363], :y)
    return ((1/4)*((y^2 - 1)^2))/10
end
    
function rho_init(t,y)
    return p_initial(y)/norminitial
end


#define the variables
#first remove any existing definitions
#unregister(model, :u)
#unregister(model, :v)
#unregister(model, :rho)


#the optimal control
@variable(model, u, Infinite(t,q), start = u_init)

#the density
@variable(model, rho>=0, Infinite(t,q), start = rho_init)

#the optimal control
@variable(model, v, Infinite(t,q), start = v_init)

#define the objective
@objective(model, Min, integral(integral((u^2)*rho,q), t)/4)

#fokker planck
@constraint(model, deriv(rho,t) - deriv(u*rho,q) - deriv(rho,q,q)== 0)
#@constraint(model, deriv(rho,t) - 10*deriv(u*rho,q) - deriv(rho,q,q)== 0)

#hjb
#@constraint(model, deriv(v,t) - 10*u*deriv(v,q) + deriv(v,q,q) + 2.5*(u^2) == 0)
@constraint(model, deriv(v,t) - u*deriv(v,q) + deriv(v,q,q) + (u^2)/4 == 0)

#stationarity condition
@constraint(model, 2*u == deriv(v,q))

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
vvals = reshape(value(v),(num_supports_t,num_supports_q));
uvals = reshape(value(u),(num_supports_t,num_supports_q));

qax = vec(unique(qgrid))

##Distribution
p2 = plot(layout = 12)

#evolve
for t in 1:11

    curr_time = plot_times[t]
    
    filtered_df = filter(row -> row.t == curr_time, df_new)
    plot!(p2,vec(qax),(vec(rhovals[t,:])),subplot=t,label = "drift output",title ="t= $curr_time",titlefontsize = 12)
    plot!(p2,vec(filtered_df.x),filtered_df.rho,subplot=t,label = "drift output")

    
end


plot!(ylim=(0,0.6))

plot!(p2,legend = false)
plot!(p2,[],[],framestyle = :none,legend = true,labels = "IPOPT", subplot = 12,fg_legend = :false)

plot!(p2,[],[],framestyle = :none,legend = true,labels = "Caluya-Halder", subplot = 12)
plot(p2)
savefig("overdamped_pdf.png")


p2 = plot(layout = 12)

#evolve

for t in 1:11

    curr_time = plot_times[t]
    
    
    filtered_df = filter(row -> row.t == curr_time, df_new)
    plot!(p2,vec(qax),(vec(uvals[t,:])),subplot=t,label = "drift output",title ="t= $curr_time",titlefontsize = 12)


    fun_temp = linear_interpolation(vec(filtered_df.x),vec(filtered_df.sigma),extrapolation_bc = Line())
    du = filtering(ForwardDiff.derivative.(Ref(fun_temp), vec(filtered_df.x)),60) #this is \partial_q V_t
    
    plot!(p2,vec(filtered_df.x),vec(du),subplot=t,label = "drift output")
    
end


plot!(p2,xlim =(-4,4))
plot!(p2,legend = false)
plot!(p2,[],[],framestyle = :none,legend = true,labels = "IPOPT", subplot = 12,fg_legend = :false)

plot!(p2,[],[],framestyle = :none,legend = true,labels = "Caluya-Halder", subplot = 12)


plot(p2)

savefig("overdamped_drift.png")


p2 = plot(layout = 12)

#evolve

for t in 1:11

    curr_time = plot_times[t]
    
    
    filtered_df = filter(row -> row.t == curr_time, df_new)
    plot!(p2,vec(qax),(vec(vvals[t,:])),subplot=t,label = "drift output",title ="t= $curr_time",titlefontsize = 12)


    #du = filtering(ForwardDiff.derivative.(Ref(fun_temp), vec(filtered_df.x)),60) #this is \partial_q V_t
    
    plot!(p2,vec(filtered_df.x),vec(filtered_df.sigma),subplot=t,label = "drift output")
    #plot!(p2, vec(qax), 0.25*(((vec(qax))).^2 .-1).^2,subplot=t)
    
end


plot!(p2,xlim =(-4,4))
plot!(p2,legend = false)
plot!(p2,[],[],framestyle = :none,legend = true,labels = "IPOPT", subplot = 12,fg_legend = :false)

plot!(p2,[],[],framestyle = :none,legend = true,labels = "Caluya-Halder", subplot = 12)


plot(p2)

savefig("overdamped_value.png")
