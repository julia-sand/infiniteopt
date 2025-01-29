using Trapz
using ArgParse

#=
This file sets up boundary conditions and parameters shared by all KL models
Arguments can be customised from the command line (parsed here)
=#


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--tf"
            help = "final time of the overdamped problem"
            arg_type = Float64
            default = 0.2
        "--epsilon"
            help = "dimensionless parameter"
            arg_type = Float64
            default = 0.2
        "--tsteps"
            help = "number of time coordinates in discretisation"
            arg_type = Int
            default = 11
        "--qsteps"
            help = "number of position coordinates in discretisation"
            arg_type = Int
            default = 10
        "--psteps"
            help = "number of momentum coordinates in discretisation"
            arg_type = Int
            default = 10
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

##LATTICE PARAMETERS!!!
T = parsed_args["tf"]
num_supports_t = parsed_args["tsteps"]
num_supports_q = parsed_args["qsteps"]
num_supports_p = parsed_args["psteps"]
epsilon = parsed_args["epsilon"]


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


global function underdamped_p_initial(p,q)

    return exp.(-((q.-1).^4)/4) .* exp.(-(p.^2)/2) /norminitial
end;


global function underdamped_p_final(p,q)

    return exp.(-(((q.^2).-1).^2)/4) .* exp.(-(p.^2)/2) /normfinal
end;
