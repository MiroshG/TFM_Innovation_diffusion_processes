using Plots, Random, DelimitedFiles, BenchmarkTools, StatsBase, Distributions, Graphs, LinearAlgebra, Profile, SpecialFunctions, DataStructures

# function to compute the time until the density of adopters reaches a certain value a=rho
function time_respect_parameters(beta, gamma,rho)
    t=1/(beta+gamma)*log((1+beta*rho/gamma)/(1-rho))
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end
#  CREATION OF NETWORKS 

function square_lattice_network(num_of_individuals,radius,L_horizontal,L_vertical)
        
    list_of_neighbors=[]
    total = 0
    for i in 1:radius-1
        total += radius-i
    end

    num_neigh=4*(radius+total)

    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for ind in 1:num_of_individuals 
        neighbors=[]
        i=trunc((ind-1)÷L_vertical+1) # row of the lattice in which the individual is (L_vertical)
        j=mod1(ind,L_horizontal)  # column in which the individual is (L_horizontal)

        # loop to obtain the neighbours 
        for r in 1:radius        # loop to obtain first neighbors, second neighbors, ...
            for col_dist in -r:r , row_dist in -r:r     # loop over the row column positions  (i,j) of the possible n neighbors
                i_neigh=mod1(i+row_dist,L_horizontal)
                j_neigh=mod1(j+col_dist,L_vertical)
                if abs(col_dist)+abs(row_dist)==r
                    push!(neighbors, Int((i_neigh-1)*L_vertical+j_neigh))
                end
            end
        end
        list_of_neighbors=vcat(list_of_neighbors,neighbors)
        dictionary_of_neighbors[ind]=neighbors
    end

    pointer_to_neighbors=[1 + num_neigh * i for i in 0:num_of_individuals]

    writedlm("networks_data//square-pointer_$num_of_individuals-$radius.txt", pointer_to_neighbors)
    writedlm("networks_data//square-list_$num_of_individuals-$radius.txt", list_of_neighbors)
    #writedlm("ratio_data//square_dict_$num_of_individuals-$radius.txt", dictionary_of_neighbors)
    writedlm("networks_data//mean_degree_square-$radius-$num_of_individuals.txt", num_neigh)

end

function random_lattice_erdos(num_of_individuals, prob)

    g=erdos_renyi(num_of_individuals, prob)
    degree_list=degree(g)
    writedlm("networks_data//mean_degree_erdos-$prob-$num_of_individuals.txt", mean(degree_list))
    list_of_neighbors=[]
    pointer_to_neighbors=[1]
    sum_neigh=1
    for i in 1:num_of_individuals
        neigh=all_neighbors(g,i)
        sum_neigh+=length(neigh)
        list_of_neighbors=vcat(list_of_neighbors,neigh)
        push!(pointer_to_neighbors, sum_neigh)
    end

    writedlm("networks_data//erdos-pointer_$num_of_individuals-$prob.txt", pointer_to_neighbors)
    writedlm("networks_data//erdos-list_$num_of_individuals-$prob.txt", list_of_neighbors)

    #=
    x=900:1100
    x_array=[i for i in x]
    lambda=num_of_individuals*prob
    poisson_dist = Poisson(lambda)
    y = pdf.(poisson_dist, x)
    histogram(degree_list, seriestype=:scatter, normalize=:pdf)
    plot!(x_array, y, linewidth=3)

    savefig("Figures\\prueba_erdos_$prob.png")
    # =#
end

function barabasi_network(num_of_individuals, degree_val)
    g=barabasi_albert(num_of_individuals, degree_val)
    degree_list=degree(g)
    writedlm("networks_data//mean_degree_RSF-$degree_val-$num_of_individuals.txt", mean(degree_list))
    list_of_neighbors=[]
    pointer_to_neighbors=[1]
    sum_neigh=1
    for i in 1:num_of_individuals
        neigh=all_neighbors(g,i)
        sum_neigh+=length(neigh)
        list_of_neighbors=vcat(list_of_neighbors,neigh)
        push!(pointer_to_neighbors, sum_neigh)
    end

    writedlm("networks_data//scale_free-pointer_$num_of_individuals-$degree_val.txt", pointer_to_neighbors)
    writedlm("networks_data//scale_free-list_$num_of_individuals-$degree_val.txt", list_of_neighbors)

    #plot(keys(degree_list), values(degree_list), seriestype=:scatter, xscale=:log10, yscale=:log10)
    #savefig("prueba.png")

    estimate_gamma(degree_list)

end

function estimate_gamma_for_scale_free(degree_list)
    degree_counts = countmap(degree_list)
    
    k_vals = collect(keys(degree_counts))  
    pk_vals = collect(values(degree_counts))  
    
    pk_vals = pk_vals ./ sum(pk_vals)
    
    k_vals = k_vals[k_vals .> 1]
    pk_vals = pk_vals[1:length(k_vals)]
    
    log_k = log.(k_vals)
    log_pk = log.(pk_vals)
    
    A = hcat(ones(length(log_k)), log_k)
    coeffs = A \ log_pk  # 
    
    gamma = -coeffs[2]
    
    println("Estimated gamma: $gamma")

end

function z_network(num_of_individuals, degree_val)
    g=random_regular_graph(num_of_individuals, degree_val)
    degree_list=degree(g)
    writedlm("networks_data//mean_degree_z-$degree_val-$num_of_individuals.txt", mean(degree_list))

    list_of_neighbors=[]
    pointer_to_neighbors=[1]
    sum_neigh=1
    for i in 1:num_of_individuals
        neigh=all_neighbors(g,i)
        sum_neigh+=length(neigh)
        list_of_neighbors=vcat(list_of_neighbors,neigh)
        push!(pointer_to_neighbors, sum_neigh)
    end

    writedlm("networks_data//z-pointer_$num_of_individuals-$degree_val.txt", pointer_to_neighbors)
    writedlm("networks_data//z-list_$num_of_individuals-$degree_val.txt", list_of_neighbors)

    #plot(keys(degree_list), values(degree_list), seriestype=:scatter, xscale=:log10, yscale=:log10)
    #savefig("prueba_z.png")

end

# ALL TO ALL 

#Discrete time algorithm for an all to all graph
function s_ratio_evolution_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,cycle)
    
    non_adopters_array=[i for i in 1:num_of_individuals]

    penetration_ratio_matrix=zeros(time+1,0)
    density_links_matrix=zeros(time+1,0)
    total_links=num_of_individuals*(num_of_individuals-1)/2

    # dictionary to store the new adopters at given times
    time_g=time/num_of_individuals
    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time_g/K:time_g
        dictionary_of_intervals[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end

    #######################################################

    montecarlo_step=1/(num_of_individuals*(gamma+beta))
    # α=γ/(γ+β)
    alpha=gamma/(gamma+beta)

    for _ in 1:cycle

        status_list=(zeros(num_of_individuals))
        penetration_ratio_list=[0.0] #sum(status_list)/num_of_individuals] #penetration_ratio sum 
        density_links_list=[0.0] 
        penetration_ratio=0.0
        density_links=0.0
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0

        for t in 0:time-1 #&& length(non_adopters_array)>0        # one-step process

            rand_elem=rand(non_adopters_array) # selection of an idividual at random
            
            if rand()<(alpha + (1-alpha)*penetration_ratio)*(1-status_list[rand_elem])
                density_links+=((num_of_individuals-1)-2*penetration_ratio*num_of_individuals)/total_links
                penetration_ratio+=1/num_of_individuals

                status_list[rand_elem]=1.0
            end
            
            push!(penetration_ratio_list, penetration_ratio)
            push!(density_links_list,density_links)

            if t*montecarlo_step>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters) # number of adopters in the interval
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end
        end

        penetration_ratio_matrix=hcat(penetration_ratio_matrix, penetration_ratio_list)
        density_links_matrix=hcat(density_links_matrix,density_links_list)
    end

    time_list=[i for i in 0:time]
    time_list=time_list.*montecarlo_step

    
    penetration_ratio_mean=[mean(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    penetration_ratio_std=[std(penetration_ratio_matrix[i,:]) for i in 1:time+1]

    density_links_mean=[mean(density_links_matrix[i,:]) for i in 1:time+1]
    density_links_std=[std(density_links_matrix[i,:]) for i in 1:time+1]

    new_adopters_mean=[]
    for i in t_inter_arr
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end
    
    results=hcat(time_list, penetration_ratio_mean)

    writedlm("ratio_data//s_ratio_evolution_Sc_n_-All_to_all-$Num-$beta-$gamma--.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma--.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)

    writedlm("ratio_data//s_ratio_evolution_Sc_n_-std-All_to_all-$Num-$beta-$gamma--.txt", result_3)

    result_4=hcat(time_list,density_links_mean)
    result_5=hcat(time_list,density_links_std)

    writedlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-All_to_all-$Num-$beta-$gamma--.txt", result_4)
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-All_to_all-$Num-$beta-$gamma--.txt", result_5)

end

#Residence time algorithm for an all to all graph
function g_ratio_evolution_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,cycle)
    
    total_links=num_of_individuals*(num_of_individuals-1)/2

    ################################## PARAMETERS ################################

    penetration=1/num_of_individuals
    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    dictionary_of_density_links=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time+1
        dictionary_of_intervals[i]=[]
        dictionary_of_density_links[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################

    for _ in 1:cycle
        
        t=0
        penetration_ratio=0
        num_individuals_array=[i for i in 1:num_of_individuals]
        status_list=(zeros(num_of_individuals))
        density_links=0.0
        i=0

        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(num_individuals_array)

            i+=1

            local_rate=gamma+beta*sum(status_list)/(length(status_list)-1)
            W_ij=local_rate*length(num_individuals_array)

            v=rand()

            wij=0
            value=0
            adopter=0
            while v*W_ij>wij
                value+=1
                wij+=local_rate
                adopter=num_individuals_array[value]
            end

            density_links+=((num_of_individuals-1)-2*penetration_ratio*num_of_individuals)/total_links
            penetration_ratio += penetration

            t_n=-log(rand())/W_ij
            t +=t_n

            status_list[adopter]=1.0

            num_individuals_array=filter(!=(adopter),num_individuals_array)

            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                push!(dictionary_of_density_links[t_inter], density_links)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end

        end
    end

    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    density_links_mean=[]
    density_links_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        push!(density_links_mean, mean(dictionary_of_density_links[i]))
        push!(density_links_std, std(dictionary_of_density_links[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)
    
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-All_to_all-$Num-$beta-$gamma--.txt", results)

    
    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma--.txt", result_2)
    

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-std-All_to_all-$Num-$beta-$gamma--.txt", result_3)

    result_4=hcat(t_inter_arr,density_links_mean)
    result_5=hcat(t_inter_arr,density_links_std)

    writedlm("ratio_data//g_ratio_evolution_Sc_n_-density_links-All_to_all-$Num-$beta-$gamma--.txt", result_4)
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-std_density_links-All_to_all-$Num-$beta-$gamma--.txt", result_5)

end

# NETWORKS 

#Discrete time algorithm for a non all to all graph
function s_ratio_evolution_Sc_n_(type,time,num_of_individuals,beta,gamma,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])

    total_links=length(list_of_neighbors)/2

    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
    end
    

    penetration_ratio_matrix=zeros(time+1,0)
    density_of_active_links_matrix=zeros(time+1,0)

    non_adopters_array=[i for i in 1:num_of_individuals]


    time_g=time/num_of_individuals
    K=100

    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time_g/K:time_g
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end
    
    montecarlo_step=1/(num_of_individuals*(gamma+beta))
    # α=γ/(γ+β)
    alpha=gamma/(gamma+beta)

    for _ in 1:cycle

        status_list=(zeros(num_of_individuals))
        penetration_ratio_list=[0.0] #penetration_ratio sum 
        density_of_active_links_list=[0.0]
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0
        density_links=0.0
        penetration=0.0

        for t in 0:time-1 # one-step process

            rand_elem=rand(non_adopters_array) # selection of an idividual at random

            # selection of all neighbors
            
            status_of_neighbors=[status_list[i] for i in dictionary_of_neighbors[rand_elem]]

            #         |[α    +   (1-α)  *                     n_v/k]                          |         (N-n)           |                          
            if rand()<(alpha + (1-alpha)*sum(status_of_neighbors)/length(status_of_neighbors))*(1-status_list[rand_elem])
                status_list[rand_elem]=1.0
                penetration+=1/num_of_individuals
                density_links+=(length(status_of_neighbors)-2*sum(status_of_neighbors))/total_links
            end

            push!(penetration_ratio_list, penetration)
            push!(density_of_active_links_list, density_links)

            if t*montecarlo_step>=t_inter
                push!(dictionary_of_new_adopters[t_inter], sum(status_list)/num_of_individuals-old_adopters)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=sum(status_list)/num_of_individuals
            end
        end

        penetration_ratio_matrix=hcat(penetration_ratio_matrix, penetration_ratio_list)
        density_of_active_links_matrix=hcat(density_of_active_links_matrix, density_of_active_links_list)

    end

    time_list=[i for i in 0:time]
    time_list=time_list.*montecarlo_step

    mean_density_of_active_links=[mean(density_of_active_links_matrix[i,:]) for i in 1:time+1]
    std_density_of_active_links=[std(density_of_active_links_matrix[i,:]) for i in 1:time+1]

    penetration_ratio_mean=[mean(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    penetration_ratio_std=[std(penetration_ratio_matrix[i,:]) for i in 1:time+1]


    results=hcat(time_list, penetration_ratio_mean)

    new_adopters_mean=[]
    for i in t_inter_arr
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end


    writedlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$Num-$beta-$gamma-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$degree_val.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt", result_3)

    result_density=hcat(time_list, mean_density_of_active_links)
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$Num-$beta-$gamma-$degree_val.txt", result_density)

    result_std_density=hcat(time_list, std_density_of_active_links)
    writedlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$degree_val.txt", result_std_density)
    
end

#Residence time algorithm for a non all to all graph
function g_ratio_evolution_Sc_n_(type,time,num_of_individuals,beta,gamma,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])
    
    total_links=length(list_of_neighbors)/2
    
    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
    end

    ################################## PARAMETERS ################################

    penetration=1/num_of_individuals
    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    dictionary_of_density_links=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time
        dictionary_of_intervals[i]=[]
        dictionary_of_density_links[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end



    ################################## MAIN LOOP ################################
    
    for _ in 1:cycle

        num_individuals_array=[i for i in 1:num_of_individuals]
        t=0.0
        penetration_ratio=0.0
        status_list=zeros(num_of_individuals)
        density_links=0.0

        w_array=ones(num_of_individuals).*gamma
        W=gamma*num_of_individuals
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0

        while t<time && !isempty(num_individuals_array)

            t_n=-log(rand())/W
            t +=t_n
            r=rand()*W
            a=0
            adopter=num_individuals_array[1]

            for n in num_individuals_array

                a+=w_array[n]

                if r<a
                    adopter=n
                    break
                end
            end

            W-=w_array[adopter]

            penetration_ratio += penetration
            status_list[adopter]=1.0

            neighbors=dictionary_of_neighbors[adopter]
            
            status_of_neighbors=[status_list[i] for i in dictionary_of_neighbors[adopter]]
            density_links+=(length(status_of_neighbors)-2*sum(status_of_neighbors))/total_links

            for neighbor in neighbors # loop over the neighbors of the adopter

                beta_neigh=beta/length(dictionary_of_neighbors[neighbor])
                s_j=status_list[neighbor]
                w_array[neighbor]+= beta_neigh*(1.0-s_j) # if it has the technology it doesn't count
                W += beta_neigh*(1.0-s_j)
            end

            num_individuals_array=filter(!=(adopter), num_individuals_array)

            #push!(time_dict[i], t)
            #push!(penetration_ratio_dict[i], penetration_ratio)
            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                push!(dictionary_of_density_links[t_inter], density_links)
                k_inter+=1
                if k_inter<=length(t_inter_arr)
                    t_inter=t_inter_arr[k_inter]
                end
                old_adopters=penetration_ratio
            end
        end

    end

    
    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    density_links_mean=[]
    density_links_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        push!(density_links_mean, mean(dictionary_of_density_links[i]))
        push!(density_links_std, std(dictionary_of_density_links[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)

    writedlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$Num-$beta-$gamma-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt", result_3)

    result_density_links=hcat(t_inter_arr, density_links_mean)
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-density_links-$type-$Num-$beta-$gamma-$degree_val.txt", result_density_links)
    
    result_std_density=hcat(t_inter_arr, density_links_std)
    writedlm("ratio_data//g_ratio_evolution_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$degree_val.txt", result_std_density)
    
end

#First time reaction algorithm for a non all to all graph
function f_ratio_evolution_Sc_n_(type,time,num_of_individuals,beta,gamma,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])

    #total_links=length(list_of_neighbors)/2

    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
    end

    ################################## PARAMETERS ################################

    penetration=1/num_of_individuals
    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    #dictionary_of_density_links=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time
        dictionary_of_intervals[i]=[]
        #dictionary_of_density_links[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################
    
    for _ in 1:cycle

        num_individuals_array=[i for i in 1:num_of_individuals]
        t=0.0
        penetration_ratio=0.0
        status_list=zeros(num_of_individuals)
        #density_links=0.0

        w_array=ones(num_of_individuals).*gamma # initial condition
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(num_individuals_array)
            
            t_n=10^8
            adopter=0

            for n in num_individuals_array

                local_rate=w_array[n]
                t_u=-log(rand())/local_rate
                if t_u<t_n
                    t_n=t_u
                    adopter=n
                end
            end

            penetration_ratio += 1/num_of_individuals
            t +=t_n
            status_list[adopter]=1.0

            neighbors=[i for i in dictionary_of_neighbors[adopter]]
            
            #status_of_neighbors=[status_list[i] for i in dictionary_of_neighbors[adopter]]
            #density_links+=(length(status_of_neighbors)-2*sum(status_of_neighbors))/total_links

            for neighbor in neighbors # loop over the neighbors of the adopter
                if status_list[neighbor]==0

                    status_of_second_neighbors=[status_list[i] for i in dictionary_of_neighbors[neighbor]]
                    w_array[neighbor]= gamma + beta*(sum(status_of_second_neighbors)/length(status_of_second_neighbors))
                end
            end

            W=sum(w_array)

            num_individuals_array=filter(!=(adopter), num_individuals_array)

            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                #push!(dictionary_of_density_links[t_inter], density_links)
                k_inter+=1
                if k_inter<=length(t_inter_arr)
                    t_inter=t_inter_arr[k_inter]
                end
                old_adopters=penetration_ratio
            end
        end

    end

    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    #density_links_mean=[]
    #density_links_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        #push!(density_links_mean, mean(dictionary_of_density_links[i]))
        #push!(density_links_std, std(dictionary_of_density_links[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)

    writedlm("qBass_data//f_ratio_evolution__Sc_n_-$type-$Num-$beta-$gamma-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("qBass_data//f_ratio_evolution_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("qBass_data//f_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt", result_3)

    #result_density_links=hcat(t_inter_arr, density_links_mean)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_density_links)
    
    #result_std_density=hcat(t_inter_arr, density_links_std)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_std_density)

end

##############################



gamma=0.0004
beta=1.088
cycles=10^3
Num=10^4
time_g=time_respect_parameters(beta, gamma, 0.999)
time=time_g*Num*(gamma+beta)

#= TYPES OF NETWORKS and the parameters selected:

-square || 1
-erdos  ||0.1
-scale_free ||5
-z  ||5; 4
-All_to_all || "-"

=#

#=
square_lattice_network(Num,1,100,100)
random_lattice_erdos(Num,0.1)
barabasi_network(Num,5)
z_network(Num,5)
z_network(Num,4)
# =#

#=
println("Discrete time algorithm")
@btime s_ratio_evolution_Sc_n_("square",Int(time),Num,beta,gamma,cycles, 1)
println("Residence time algorithm")
@btime g_ratio_evolution_Sc_n_("square",time_g,Num,beta,gamma,cycles, 1)
println("First reactio algorithm")
@btime f_ratio_evolution_Sc_n_("square",time_g,Num,beta,gamma,cycles, 1)
# =#

#s_ratio_evolution_Sc_n_All_to_all(Int(round(time)),Num,beta,gamma,cycles)
#g_ratio_evolution_Sc_n_All_to_all(time_g,Num,beta,gamma,cycles)

#s_ratio_evolution_Sc_n_("square",Int(round(time)),Num,beta,gamma,cycles, 1)
#g_ratio_evolution_Sc_n_("square",time_g,Num,beta,gamma,cycles, 1)

#s_ratio_evolution_Sc_n_("scale_free",Int(round(time)),Num,beta,gamma,cycles, 5)
#g_ratio_evolution_Sc_n_("scale_free",time_g,Num,beta,gamma,cycles, 5)

#s_ratio_evolution_Sc_n_("z",Int(round(time)),Num,beta,gamma,cycles, 5)
#g_ratio_evolution_Sc_n_("z",time_g,Num,beta,gamma,cycles, 5)

#s_ratio_evolution_Sc_n_("z",Int(round(time)),Num,beta,gamma,cycles, 4)
#g_ratio_evolution_Sc_n_("z",time_g,Num,beta,gamma,cycles, 4)

#s_ratio_evolution_Sc_n_("erdos",Int(round(time)),Num,beta,gamma,cycles, 0.1)
#g_ratio_evolution_Sc_n_("erdos",time_g,Num,beta,gamma,cycles, 0.1)

################################### PLOTS ####################################

function plotting_adopters_gaussian_std(type,N,beta,gamma,degree_val,time,cycles)

    result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
    std_1=result_1[:,2]./sqrt(2*cycles)

    x_data=result_1[:,1]
    y_data=result_1[:,2]
    y_std=std_1
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    gaussian_result=readdlm("ratio_data//Gaussian_approx_std_$beta-$gamma.txt")

    plot(result_1[:,1],result_1[:,2], label="Discrete time algorithm", linewidth=4, linecolor=:gray, grid=false, xlabel="Time (years)", ylabel="Standard deviation")
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:gray)
    plot!(gaussian_result[:,1],gaussian_result[:,2], label="Gaussian approximation", linestyle=:dot, linewidth=4, linecolor=:black)
    xlims!(0,time)
    savefig("Figures//Gaussian//std-$type-$N-$beta-$gamma-$degree_val.png")
    
end

function plotting_all_together_with_SAME(type,beta,gamma,N,degree_val,time)
    if type=="square"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_square-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Pair_approx-square-$N-$beta-$gamma-$mu.txt")
    elseif type=="z"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_z-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Pair_approx-z-$N-$beta-$gamma-$mu.txt")
    end
        
    AME_result=readdlm("ratio_data//AME-$type-$N-$beta-$gamma-$degree_val.txt")
    gaussian_result=readdlm("ratio_data//Gaussian_approx_mean_$beta-$gamma.txt")
    SAME_result=readdlm("ratio_data/SAME-adopters-$type-$N-$beta-$gamma-$degree_val.txt")
    mf_result=readdlm("ratio_data//MF_approx_mean_$beta-$gamma.txt")

    x_data=result_1[:,1]
    y_data=result_1[:,2]
    y_std=std_1[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_1[:,1],result_1[:,2], label="Discrete time algorithm", linewidth=4, linecolor=:gray, grid=false, xlabel="Time (years)", ylabel="Density of adopters")
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:gray)
    #plot!(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2)
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(gaussian_result[:,1],gaussian_result[:,2], label="Gaussian approximation", linestyle=:dot, linewidth=4, linecolor=:black)
    plot!(pair_approx[:,1],pair_approx[:,2], label="Pair approximation", linewidth=4, linecolor=:red)
    plot!(AME_result[:,1],AME_result[:,2], label="AME", linestyle=:dot, linewidth=4, linecolor=:blue)
    plot!(SAME_result[:,1],SAME_result[:,2], label="SAME", linestyle=:dot, linewidth=4, linecolor=:orange)
    xlims!(0,time)

    savefig("Figures//Bass_model//Comparison_every_approximation-$type-$N-$beta-$gamma-$degree_val.png")
end

function plotting_all_together_without_SAME(type,beta,gamma,N,degree_val,time)
    if type=="erdos"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_erdos-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Pair_approx-erdos-$N-$beta-$gamma-$mu.txt")
    elseif type=="scale_free"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_RSF-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Pair_approx-RSF-$N-$beta-$gamma-$mu.txt")
    elseif type=="All_to_all"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu=9999.0
        pair_approx=readdlm("ratio_data//Pair_approx-All_to_all-$N-$beta-$gamma-$mu.txt")
    end
        
    AME_result=readdlm("ratio_data//AME-$type-$N-$beta-$gamma-$degree_val.txt")
    gaussian_result=readdlm("ratio_data//Gaussian_approx_mean_$beta-$gamma.txt")
    mf_result=readdlm("ratio_data//MF_approx_mean_$beta-$gamma.txt")

    x_data=result_1[:,1]
    y_data=result_1[:,2]
    y_std=std_1[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_1[:,1],result_1[:,2], label="Discrete time algorithm", linewidth=4, linecolor=:gray, grid=false, xlabel="Time (years)", ylabel="Density of adopters")
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:gray)
    #plot!(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2)
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(gaussian_result[:,1],gaussian_result[:,2], label="Gaussian approximation", linestyle=:dot, linewidth=4, linecolor=:black)
    plot!(pair_approx[:,1],pair_approx[:,2], label="Pair approximation", linewidth=4, linecolor=:red)
    plot!(AME_result[:,1],AME_result[:,2], label="AME", linestyle=:dot, linewidth=4, linecolor=:blue)
    #plot!(SAME_result[:,1],SAME_result[:,2], label="SAME", linestyle=:dot, linewidth=4, linecolor=:orange)
    xlims!(0,time)

    savefig("Figures//Bass_model//Comparison_every_approximation-$type-$N-$beta-$gamma-$degree_val.png")
end

function plotting_all_together_active_links_with_SAME(type,beta,gamma,N,degree_val,time)
    if type=="square"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_square-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Density_approx-square-$N-$beta-$gamma-$mu.txt")
    elseif type=="z"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_z-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Density_approx-z-$N-$beta-$gamma-$mu.txt")
    end
        
    SAME_result=readdlm("ratio_data//SAME-links-$type-$beta-$gamma-$degree_val.txt")
    AME_result=readdlm("ratio_data//AME-density_links-$type-$N-$beta-$gamma-$degree_val.txt")

    x_data=result_1[:,1]
    y_data=result_1[:,2]
    y_std=std_1[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_1[:,1],result_1[:,2], label="Discrete time algorithm", linewidth=4, linecolor=:gray, grid=false, xlabel="Time (years)", ylabel="Density of active links",legend=:topleft)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:gray)
    #plot!(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2)
    plot!(pair_approx[:,1],pair_approx[:,2], label="Pair approximation", linewidth=4, linecolor=:red)
    plot!(AME_result[:,1],AME_result[:,2], label="AME", linestyle=:dot, linewidth=4, linecolor=:blue)
    plot!(SAME_result[:,1],SAME_result[:,2], label="SAME", linestyle=:dot, linewidth=4, linecolor=:orange)
    xlims!(0,time)

    savefig("Figures//Bass_model//Comparison_every_approximation_density_links-$type-$N-$beta-$gamma-$degree_val.png")
end

function plotting_all_together_active_links_without_SAME(type,beta,gamma,N,degree_val,time)
    if type=="erdos"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_erdos-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Density_approx-erdos-$N-$beta-$gamma-$mu.txt")
    elseif type=="scale_free"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu_arr=readdlm("ratio_data//mean_degree_RSF-$degree_val.txt")[1:1]
        mu=mu_arr[1]
        pair_approx=readdlm("ratio_data//Density_approx-RSF-$N-$beta-$gamma-$mu.txt")
    elseif type=="All_to_all"
        result_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #result_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-$type-$N-$beta-$gamma-$degree_val.txt")
        std_1=readdlm("ratio_data//s_ratio_evolution_Sc_n_-std_density_links-$type-$N-$beta-$gamma-$degree_val.txt")
        #std_2=readdlm("ratio_data//g_ratio_evolution_Sc_n_-std-$type-$Num-$beta-$gamma-$degree_val.txt")
        mu=9999.0
        pair_approx=readdlm("ratio_data//Density_approx-All_to_all-$N-$beta-$gamma-$mu.txt")
    end
        
    AME_result=readdlm("ratio_data//AME-density_links-$type-$N-$beta-$gamma-$degree_val.txt")

    x_data=result_1[:,1]
    y_data=result_1[:,2]
    y_std=std_1[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_1[:,1],result_1[:,2], label="Discrete time algorithm", linewidth=4, linecolor=:gray, grid=false, xlabel="Time (years)", ylabel="Density of active links",legend=:topleft)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:gray)
    #plot!(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2)
    plot!(pair_approx[:,1],pair_approx[:,2], label="Pair approximation", linewidth=4, linecolor=:red)
    plot!(AME_result[:,1],AME_result[:,2], label="AME", linestyle=:dot, linewidth=4, linecolor=:blue)
    xlims!(0,time)

    savefig("Figures//Bass_model//Comparison_every_approximation_density_links-$type-$N-$beta-$gamma-$degree_val.png")
end
###############################################################################

#=
plotting_adopters_gaussian_std("All_to_all",Num,beta,gamma,"-",time_g,cycles)
plotting_adopters_gaussian_std("square",Num,beta,gamma,1,time_g,cycles)
plotting_adopters_gaussian_std("scale_free",Num,beta,gamma,5,time_g,cycles)
plotting_adopters_gaussian_std("z",Num,beta,gamma,5,time_g,cycles)
plotting_adopters_gaussian_std("z",Num,beta,gamma,4,time_g,cycles)
plotting_adopters_gaussian_std("erdos",Num,beta,gamma,0.1,time_g,cycles)
# =#

#=
plotting_all_together_with_SAME("square",beta,gamma,Num,1,time_g)
plotting_all_together_with_SAME("z",beta,gamma,Num,4,time_g)
plotting_all_together_with_SAME("z",beta,gamma,Num,5,time_g)
plotting_all_together_without_SAME("scale_free",beta,gamma,Num,5,time_g)
plotting_all_together_without_SAME("erdos",beta,gamma,Num,0.1,time_g)
plotting_all_together_without_SAME("All_to_all",beta,gamma,Num,"-",time_g)
# =#

#=
plotting_all_together_active_links_with_SAME("square",beta,gamma,Num,1,time_g)
plotting_all_together_active_links_with_SAME("z",beta,gamma,Num,4,time_g)
plotting_all_together_active_links_with_SAME("z",beta,gamma,Num,5,time_g)
plotting_all_together_active_links_without_SAME("scale_free",beta,gamma,Num,5,time_g)
plotting_all_together_active_links_without_SAME("erdos",beta,gamma,Num,0.1,time_g)
#plotting_all_together_active_links_without_SAME("All_to_all",beta,gamma,Num,"-",time_g)
# =#
