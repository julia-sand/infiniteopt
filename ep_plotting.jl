using Plots;
using CSV;
using DataFrames;
#using ImageFiltering;
#using Interpolations;
#using ForwardDiff;

T = 0.2

##filtering functions
function filtering(input_array,filter_delta)
    
    #= this function implements a simple box filter of size filter_delta and padds missing reults
    
    inputs: 
    input_array - array to be smooths
    filter_delta - size of filter, number of nearest neighbours to consider
    
    output:
    filtered array with borders filled with constant value. 
    =#
    
    return imfilter(vec(input_array),reflect((1/filter_delta).*centered(ones(filter_delta))),Pad(:replicate))
    
end;

#get both dataframes

#get Caluya-Halder output
#df_calhal = CSV.read("infiniteopt/results_kl_land_new1.csv",DataFrame,header=true);
df_ipopt = CSV.read("infiniteopt/ep/ipopt_overdampedep_v1.csv",DataFrame,header=true);

plot_times = range(0, T, 11)


p2 = plot(layout = 12)

##density plots
for t in 1:11

    curr_time = plot_times[t]
    println(curr_time)    
    
    #filtered_df = filter(row -> row.t == curr_time, df_calhal)
    filtered_df_ipopt = filter(row -> row.t == curr_time, df_ipopt)


    #plot!(p2,vec(filtered_df.x),vec(filtered_df.rho),subplot=t,label = "Caluya-Halder")
    plot!(p2,vec(filtered_df_ipopt.x),vec(filtered_df_ipopt.rho),subplot=t,label = "InfiniteOpt",title ="t= $curr_time",titlefontsize = 12)
    
end


plot!(p2,xlim =(-4,4))
plot!(p2,ylim =(-0.05,0.5))
plot!(p2,legend = false)

#plot!(p2,[],[],framestyle = :none,legend = true,labels = "Caluya-Halder", subplot = 12,fg_legend = :false)
plot!(p2,[],[],framestyle = :none,legend = true,labels = "IPOPT", subplot = 12)


#plot(p2)

savefig("ep_overdamped_distributions.png")

