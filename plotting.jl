using Plots;

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
