using Plots, Random, DelimitedFiles, BenchmarkTools, StatsBase, Distributions, Graphs, LinearAlgebra, Profile, SpecialFunctions, DataStructures

function time_respect_parameters(beta, gamma,q,rho)

    mf_result=readdlm("ratio_data//MF_qBass_approx_mean_$beta-$gamma-$q.txt")
    position=findfirst(x->x>rho, mf_result[:,2])
    t=mf_result[:,1][position]
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end

# ALL TO ALL 

#Discrete time algorithm for an all to all graph
function s_ratio_evolution_qBass_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,q,cycle)
    
    non_adopters_array=[i for i in 1:num_of_individuals]

    penetration_ratio_matrix=zeros(time+1,0)

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
        penetration_ratio=0.0
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0

        for t in 0:time-1 #&& length(non_adopters_array)>0        # one-step process

            rand_elem=rand(non_adopters_array) # selection of an idividual at random
            
            if rand()<(alpha + (1-alpha)*penetration_ratio^q)*(1-status_list[rand_elem])
                status_list[rand_elem]=1.0
                penetration_ratio+=1/num_of_individuals
            end
            
            push!(penetration_ratio_list, penetration_ratio)

            if t*montecarlo_step>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters) # number of adopters in the interval
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end
        end

        penetration_ratio_matrix=hcat(penetration_ratio_matrix, penetration_ratio_list)
    end

    time_list=[i for i in 0:time]
    time_list=time_list.*montecarlo_step

    
    penetration_ratio_mean=[mean(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    penetration_ratio_std=[std(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    
    new_adopters_mean=[]
    for i in t_inter_arr
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end
    
    results=hcat(time_list, penetration_ratio_mean)

    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-All_to_all-$Num-$beta-$gamma-$q.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$q.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)

    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$q.txt", result_3)

end

#Residence time algorithm for an all to all graph
function g_ratio_evolution_qBass_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,q,cycle)
    

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time+1
        dictionary_of_intervals[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################

    for _ in 1:cycle
        
        t=0
        penetration_ratio=0
        num_individuals_array=[i for i in 1:num_of_individuals]
        status_list=(zeros(num_of_individuals))

        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(num_individuals_array)


            local_rate=gamma+beta*(sum(status_list)/(length(status_list)-1))^q
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

            penetration_ratio += 1/num_of_individuals
            t_n=-log(rand())/W_ij
            t +=t_n

            status_list[adopter]=1.0

            num_individuals_array=filter(!=(adopter),num_individuals_array)

            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end

        end
    end

    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)
    
    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-All_to_all-$Num-$beta-$gamma-$q.txt", results)

    
    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$q.txt", result_2)
    

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$q.txt", result_3)

end

#First time reaction algorithm for an all to all graph
function f_ratio_evolution_qBass_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,q,cycle)

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time+1
        dictionary_of_intervals[i]=[]
        dictionary_of_new_adopters[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################

    for _ in 1:cycle

        t=0
        penetration_ratio=0
        num_individuals_array=[i for i in 1:num_of_individuals]
        #status_list=(zeros(num_of_individuals))

        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(num_individuals_array)
            

            t_n=10^8
            adopter=0

            local_rate=gamma+beta*penetration_ratio^q #(sum(status_list)/(length(status_list)-1))^q

            for n in num_individuals_array

                t_u=-log(rand())/local_rate
                if t_u<t_n
                    t_n=t_u
                    adopter=n
                end
            end
            penetration_ratio += 1/num_of_individuals
            t +=t_n
            status_list[adopter]=1.0

            num_individuals_array=filter(!=(adopter), num_individuals_array)

            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end
        end

    end

    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)
    
    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-All_to_all-$Num-$beta-$gamma-$q.txt", results)
    
    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$q.txt", result_2)
    

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$q.txt", result_3)
end

# NETWORKS 

#Discrete time algorithm for a non all to all graph
function s_ratio_evolution_qBass_Sc_n_(type,time,num_of_individuals,beta,gamma,q,cycle,degree_val)
    
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

            ratio_to_the_square=1
            for neigh in 1:q
                rand_neighbor=rand(status_of_neighbors)
                ratio_to_the_square*=rand_neighbor
            end

            #         |[α    +   (1-α)  *       n_v/k]       |         (N-n)           |                          
            if rand()<(alpha + (1-alpha)*ratio_to_the_square)*(1-status_list[rand_elem])
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


    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-$type-$Num-$beta-$gamma-$q-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)
    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-std-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_3)

    result_density=hcat(time_list, mean_density_of_active_links)
    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_density)

    result_std_density=hcat(time_list, std_density_of_active_links)
    writedlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_std_density)
    
end

#Residence time algorithm for a non all to all graph
function g_ratio_evolution_qBass_Sc_n_(type,time,num_of_individuals,beta,gamma,q,cycle,degree_val)
    
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

            penetration_ratio += penetration
            status_list[adopter]=1.0
            w_array[adopter]=0.0 # probability of adoption of an adopter is 0

            neighbors=dictionary_of_neighbors[adopter]
            
            #status_of_neighbors=[status_list[i] for i in dictionary_of_neighbors[adopter]]
            #density_links+=(length(status_of_neighbors)-2*sum(status_of_neighbors))/total_links

            for neighbor in neighbors # loop over the neighbors of the adopter
                if status_list[neighbor]==0

                    status_of_second_neighbors=[status_list[i] for i in dictionary_of_neighbors[neighbor]]
                    w_array[neighbor]= gamma + (sum(status_of_second_neighbors)/length(status_of_second_neighbors))^q 
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

    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-$type-$Num-$beta-$gamma-$q-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_3)

    #result_density_links=hcat(t_inter_arr, density_links_mean)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_density_links)
    
    #result_std_density=hcat(t_inter_arr, density_links_std)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_std_density)
    
end

#First time reaction algorithm for a non all to all graph
function f_ratio_evolution_qBass_Sc_n_(type,time,num_of_individuals,beta,gamma,q,cycle,degree_val)
    
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
                    w_array[neighbor]= gamma + beta*(sum(status_of_second_neighbors)/length(status_of_second_neighbors))^q 
                end
            end

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

    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-$type-$Num-$beta-$gamma-$q-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("qBass_data//f_ratio_evolution_qBass_Sc_n_-std-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_3)

    #result_density_links=hcat(t_inter_arr, density_links_mean)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_density_links)
    
    #result_std_density=hcat(t_inter_arr, density_links_std)
    #writedlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$q-$degree_val.txt", result_std_density)

end


##############################



gamma=0.0004 #0.0004 
beta=1.088 #1.088 
#alpha=0.1 # alpha < 1
q=2
cycles=10^2
Num=10^4
time_g=100 #time_respect_parameters(beta, gamma, q, 0.999)
time=time_g*Num*(gamma+beta)

#= TYPES OF NETWORKS:

-square || 1
-erdos  ||0.1
-scale_free ||5
-z  ||5; 4

=#

#=
println("Discrete")
@btime s_ratio_evolution_qBass_Sc_n_("square",Int(round(time)),Num,beta,gamma,q,cycles, 1)
println("Residence time")
@btime g_ratio_evolution_qBass_Sc_n_("square",time_g,Num,beta,gamma,q,cycles, 1)
println("First reaction")
@btime f_ratio_evolution_qBass_Sc_n_("square",time_g,Num,beta,gamma,q,cycles, 1)
=#

#######################################################################

function plotting_adopters_MF_qBass_All_to_all(N,beta,gamma,q,time)


    #result_1=readdlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-All_to_all-$N-$beta-$gamma-$q.txt")
    result_2=readdlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-All_to_all-$N-$beta-$gamma-$q.txt")
    std_2=readdlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std-All_to_all-$N-$beta-$gamma-$q.txt")
    

    mf_result=readdlm("qBass_data//MF_qBass_approx_mean_$beta-$gamma-$q.txt")
    
    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2, grid=false,xlabel="Time (years)",ylabel="Density of adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//qBass//Adopters_all_to_all-$N-$beta-$gamma-$q.png")

end

#######################################################################

function plotting_adopters_MF_qBass(type,N,beta,gamma,q,degree_val,time)


    #result_1=readdlm("qBass_data//s_ratio_evolution_qBass_Sc_n_-$type-$N-$beta-$gamma-$q-$degree_val.txt")
    result_2=readdlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-$type-$N-$beta-$gamma-$q-$degree_val.txt")
    std_2=readdlm("qBass_data//g_ratio_evolution_qBass_Sc_n_-std-$type-$N-$beta-$gamma-$q-$degree_val.txt")

    mf_result=readdlm("qBass_data//MF_qBass_approx_mean_$beta-$gamma-$q.txt")
        
    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2, grid=false,xlabel="Time (years)",ylabel="Density of adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//qBass//Adopters_$type-$N-$beta-$gamma-$q-$degree_val.png")
end

#######################################################################

#s_ratio_evolution_qBass_Sc_n_All_to_all(Int(round(time)),Num,beta,gamma,q,cycles)
#g_ratio_evolution_qBass_Sc_n_All_to_all(time_g,Num,beta,gamma,q,cycles)



#=
println("square")
g_ratio_evolution_qBass_Sc_n_("square",time_g,Num,beta,gamma,q,cycles, 1)

println("scale_free")
g_ratio_evolution_qBass_Sc_n_("scale_free",time_g,Num,beta,gamma,q,cycles, 5)

println("z 5")
g_ratio_evolution_qBass_Sc_n_("z",time_g,Num,beta,gamma,q,cycles, 5)

println("z 4")
g_ratio_evolution_qBass_Sc_n_("z",time_g,Num,beta,gamma,q,cycles, 4)

println("erdos")
g_ratio_evolution_qBass_Sc_n_("erdos",time_g,Num,beta,gamma,q,cycles, 0.1)
# =#

# =
plotting_adopters_MF_qBass_All_to_all(Num,beta,gamma,q,time_g)
plotting_adopters_MF_qBass("square",Num,beta,gamma,q,1,time_g)
plotting_adopters_MF_qBass("scale_free",Num,beta,gamma,q,5,time_g)
plotting_adopters_MF_qBass("z",Num,beta,gamma,q,5,time_g)
plotting_adopters_MF_qBass("z",Num,beta,gamma,q,4,time_g)
plotting_adopters_MF_qBass("erdos",Num,beta,gamma,q,0.1,time_g)
# =#