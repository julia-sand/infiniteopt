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
        "--tsteps"
            help = "number of time coordinates in discretisation"
            arg_type = Int
            default = 11
        "--qsteps"
            help = "number of space coordinates in discretisation"
            arg_type = Int
            default = 20
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()


##LATTICE PARAMETERS!!!
T = 0.2
num_supports_t = 101
num_supports_q = 1000



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

