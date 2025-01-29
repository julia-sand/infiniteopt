using Trapz;
#=
This file sets up boundary conditions and parameters shared by all KL models
=#


#####boundary conditions#######
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
T = 0.2
num_supports_t = 101
num_supports_q = 1000
