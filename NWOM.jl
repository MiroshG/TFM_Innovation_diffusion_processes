using Plots, Random, DelimitedFiles, BenchmarkTools, StatsBase, Distributions, Graphs, LinearAlgebra, Profile, SpecialFunctions, DataStructures

function time_respect_parameters_NWOM(beta, gamma,delta,alpha,rho)
    
    mf_result=readdlm("ratio_data//MF_NWOM_adopters_$beta-$gamma-$delta-$alpha.txt")

    position=findfirst(x->x>rho, mf_result[:,2])
    t=mf_result[:,1][position]
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end

# ALL TO ALL 

function s_ratio_evolution_NWOM_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,delta,alpha,cycle)
    
    non_adopters_array=[i for i in 1:num_of_individuals]

    penetration_ratio_matrix=zeros(time+1,0)
    hate_ratio_matrix=zeros(time+1,0)

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
    eta=gamma/(gamma+beta)

    for _ in 1:cycle

        status_list=(zeros(num_of_individuals))
        penetration_ratio_list=[0.0] #sum(status_list)/num_of_individuals] #penetration_ratio sum 
        hate_ratio_list=[0.0]
        penetration_ratio=0.0
        hate_ratio=0.0
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0

        for t in 0:time-1 #&& length(non_adopters_array)>0        # one-step process

            rand_elem=rand(non_adopters_array) # selection of an idividual at random
            
            if status_list[rand_elem]==0 # selecting a non adopter
                u=rand()
                if u<delta*(eta+(1-eta)*(penetration_ratio-hate_ratio))*(1-hate_ratio)^alpha #hates
                    hate_ratio+=1/num_of_individuals
                    penetration_ratio+=1/num_of_individuals # either if it hates the product or not, it have adopted
                    status_list[rand_elem]=-1.0

                elseif u<(eta+(1-eta)*(penetration_ratio-hate_ratio))*(1-hate_ratio)^alpha #happy adopter
                    penetration_ratio+=1/num_of_individuals
                    status_list[rand_elem]=1.0
                end
            end
            
            push!(penetration_ratio_list, penetration_ratio)
            push!(hate_ratio_list, hate_ratio)

            if t*montecarlo_step>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters) # number of adopters in the interval
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end
        end

        penetration_ratio_matrix=hcat(penetration_ratio_matrix, penetration_ratio_list)
        hate_ratio_matrix=hcat(hate_ratio_matrix, hate_ratio_list)
    end

    time_list=[i for i in 0:time]
    time_list=time_list.*montecarlo_step

    
    penetration_ratio_mean=[mean(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    penetration_ratio_std=[std(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    hate_ratio_mean=[mean(hate_ratio_matrix[i,:]) for i in 1:time+1]
    hate_ratio_std=[std(hate_ratio_matrix[i,:]) for i in 1:time+1]
    

    new_adopters_mean=[]
    for i in t_inter_arr
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end
    
    results=hcat(time_list, penetration_ratio_mean)

    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)

    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_3)

    result_4=hcat(time_list, hate_ratio_mean)

    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_4)

    result_5=hcat(time_list, hate_ratio_std)

    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_5)


end

function g_ratio_evolution_NWOM_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,delta,alpha,cycle)
    

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    dictionary_of_hate_intervals=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time+1
        dictionary_of_intervals[i]=[]
        dictionary_of_new_adopters[i]=[]
        dictionary_of_hate_intervals[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################

    for _ in 1:cycle
        
        non_adopters_array=[i for i in 1:num_of_individuals] # array with the non adopters
        posibilities_array=[i for i in 1:2*num_of_individuals] # possible transitions for the non adopters (2 times the number of non adopters)
        t=0.0
        penetration_ratio=0.0
        hate_ratio=0.0

        w_array_adopt=ones(num_of_individuals).*gamma*(1-delta) #positive adoption
        w_array_hates=ones(num_of_individuals).*gamma*delta     #negative adoption
        w_arrays=vcat(w_array_adopt,w_array_hates)
        W=num_of_individuals*gamma # initial chances to adopt
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(posibilities_array)


            t_n=-log(rand())/W
            t +=t_n
            r=rand()*W
            a=0
            adopter=posibilities_array[1]
            i=0

            for n in posibilities_array

                i+=1
                a+=w_arrays[n]

                if r<a
                    adopter=n
                    break
                end
            end

            if i>length(non_adopters_array) # the adopter hates the technology
                hate_ratio+=1/num_of_individuals

                adopter_2=Int(adopter-num_of_individuals) # position of the possibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of unhappy adoption
                non_adopters_array=filter(!=(adopter_2),non_adopters_array) # we erase the adopter from the non adopters list

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0

            else # the adopter is happy with the technology

                adopter_2=Int(adopter+num_of_individuals) # position of the possibility of unhappy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of unhappy adoption
                non_adopters_array=filter(!=(adopter),non_adopters_array) # we erase the adopter from the non adopters list

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0
            end

            penetration_ratio+=1/num_of_individuals

            prob_copy=(gamma+beta*(penetration_ratio-hate_ratio))*(1-hate_ratio)^alpha

            for i in non_adopters_array # change in the prob of every non adopter to positive or negative adopt
                w_arrays[i]=prob_copy*(1-delta)                     # positive adopt 
                w_arrays[Int(i+num_of_individuals)]=prob_copy*delta # negative adopt
            end
            
            W=sum(w_arrays)


            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                push!(dictionary_of_hate_intervals[t_inter], hate_ratio)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=penetration_ratio
            end

        end
    end

    penetration_mean=[]
    penetration_std=[]
    new_adopters_mean=[]
    hate_mean=[]
    hate_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        push!(hate_mean, mean(dictionary_of_hate_intervals[i]))
        push!(hate_std, std(dictionary_of_hate_intervals[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)
    
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", results)

    
    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_2)
    

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_3)

    result_4=hcat(t_inter_arr, hate_mean)
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_4)

    result_5=hcat(t_inter_arr, hate_std)
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_5)

end

function f_ratio_evolution_NWOM_Sc_n_All_to_all(time,num_of_individuals,beta,gamma,delta,alpha,cycle)

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    dictionary_of_hate_intervals=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time+1
        dictionary_of_intervals[i]=[]
        dictionary_of_new_adopters[i]=[]
        dictionary_of_hate_intervals[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################

    for _ in 1:cycle

        non_adopters_array=[i for i in 1:num_of_individuals] # array with the non adopters
        posibilities_array=[i for i in 1:2*num_of_individuals] # possible transitions for the non adopters (2 times the number of non adopters)
        t=0.0
        penetration_ratio=0.0
        hate_ratio=0.0

        w_array_adopt=ones(num_of_individuals).*gamma*(1-delta) #positive adoption
        w_array_hates=ones(num_of_individuals).*gamma*delta     #negative adoption
        w_arrays=vcat(w_array_adopt,w_array_hates)
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(num_individuals_array)

            t_n=10^8
            adopter=0

            i=0
            for n in posibilities_array

                local_rate=w_arrays[n]
                t_u=-log(rand())/local_rate
                i+=1
                if t_u<t_n
                    t_n=t_u
                    adopter=n
                end
            end

            if i>length(non_adopters_array) # the adopter hates the technology
                hate_ratio+=1/num_of_individuals

                adopter_2=Int(adopter-num_of_individuals) # position of the possibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of unhappy adoption
                non_adopters_array=filter(!=(adopter_2),non_adopters_array) # we erase the adopter from the non adopters list

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0

            else # the adopter is happy with the technology

                adopter_2=Int(adopter+num_of_individuals) # position of the possibility of unhappy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of unhappy adoption
                non_adopters_array=filter(!=(adopter),non_adopters_array) # we erase the adopter from the non adopters list

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0
            end
            penetration_ratio+=1/num_of_individuals

            prob_copy=(gamma+beta*(penetration_ratio-hate_ratio))*(1-hate_ratio)^alpha

            for i in non_adopters_array # change in the prob of every non adopter to positive or negative adopt
                w_arrays[i]=prob_copy*(1-delta)                     # positive adopt 
                w_arrays[Int(i+num_of_individuals)]=prob_copy*delta # negative adopt
            end
            
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
    hate_mean=[]
    hate_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        push!(hate_mean, mean(dictionary_of_hate_intervals[i]))
        push!(hate_std, std(dictionary_of_hate_intervals[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)
    
    writedlm("ratio_data//f_ratio_evolution_NWOM_Sc_n_-All_to_all-$Num-$beta-$gamma-$q.txt", results)
    
    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("ratio_data//f_ratio_evolution_NWOM_Sc_n_-new_adopt-All_to_all-$Num-$beta-$gamma-$q.txt", result_2)
    

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("ratio_data//f_ratio_evolution_NWOM_Sc_n_-std-All_to_all-$Num-$beta-$gamma-$q.txt", result_3)

    result_4=hcat(t_inter_arr, hate_mean)
    writedlm("ratio_data//f_ratio_evolution_NWOM_Sc_n_-hate-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_4)

    result_5=hcat(t_inter_arr, hate_std)
    writedlm("ratio_data//f_ratio_evolution_NWOM_Sc_n_-hate-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt", result_5)
end

# NETWORKS 

function s_ratio_evolution_NWOM_Sc_n_(type,time,num_of_individuals,beta,gamma,delta,alpha,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])

    total_links=length(list_of_neighbors)/2

    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
    end
    

    penetration_ratio_matrix=zeros(time+1,0)
    hate_ratio_matrix=zeros(time+1,0)
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
    eta=gamma/(gamma+beta)

    for _ in 1:cycle

        status_list=(zeros(num_of_individuals))
        penetration_ratio_list=[0.0] #penetration_ratio sum 
        hate_ratio_list=[0.0]
        density_of_active_links_list=[0.0]
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0
        density_links=0.0
        penetration=0.0
        hate_ratio=0.0

        for t in 0:time-1 # one-step process

            rand_elem=rand(non_adopters_array) # selection of an idividual at random

            #=
            counter_status=counter(status_of_neighbors)

            # number of positive and negative neighbors:
            
            positive_density=counter_status[1.0]/length(status_of_neighbors)
            negative_density=counter_status[-1.0]/length(status_of_neighbors)
            =#

            if status_list[rand_elem]==0 # selecting a non adopter
                u=rand()
                
                # selection of all neighbors
            
                status_of_neighbors=[status_list[i] for i in dictionary_of_neighbors[rand_elem]]
                
                positive_density=count(x -> x == 1.0, status_of_neighbors)/length(status_of_neighbors)  # n_+/k
                negative_density=count(x -> x == -1.0, status_of_neighbors)/length(status_of_neighbors) # n_-/k

                if u <(eta+(1-eta)*positive_density)*(1-negative_density)^alpha*delta # adopts and hates
                    status_list[rand_elem]=-1.0 # hates
                    penetration+=1/num_of_individuals
                    hate_ratio+=1/num_of_individuals
                elseif u <(eta+(1-eta)*positive_density)*(1-negative_density)^alpha   # only adopts
                    status_list[rand_elem]=1.0 # adopts
                    penetration+=1/num_of_individuals
                end
            end

            push!(penetration_ratio_list, penetration)
            push!(density_of_active_links_list, density_links)
            push!(hate_ratio_list, hate_ratio)

            if t*montecarlo_step>=t_inter
                push!(dictionary_of_new_adopters[t_inter], sum(status_list)/num_of_individuals-old_adopters)
                k_inter+=1
                t_inter=t_inter_arr[k_inter]
                old_adopters=sum(status_list)/num_of_individuals
            end
        end

        penetration_ratio_matrix=hcat(penetration_ratio_matrix, penetration_ratio_list)
        density_of_active_links_matrix=hcat(density_of_active_links_matrix, density_of_active_links_list)
        hate_ratio_matrix=hcat(hate_ratio_matrix,hate_ratio_list)

    end

    time_list=[i for i in 0:time]
    time_list=time_list.*montecarlo_step

    mean_density_of_active_links=[mean(density_of_active_links_matrix[i,:]) for i in 1:time+1]
    std_density_of_active_links=[std(density_of_active_links_matrix[i,:]) for i in 1:time+1]

    penetration_ratio_mean=[mean(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    penetration_ratio_std=[std(penetration_ratio_matrix[i,:]) for i in 1:time+1]
    hate_ratio_mean=[mean(hate_ratio_matrix[i,:]) for i in 1:time+1]
    hate_ratio_std=[std(hate_ratio_matrix[i,:]) for i in 1:time+1]


    results=hcat(time_list, penetration_ratio_mean)

    new_adopters_mean=[]
    for i in t_inter_arr
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
    end


    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time_g))
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_2)

    result_3=hcat(time_list, penetration_ratio_std)
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_3)

    result_density=hcat(time_list, mean_density_of_active_links)
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_density)

    result_std_density=hcat(time_list, std_density_of_active_links)
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_density)

    result_hate=hcat(time_list, hate_ratio_mean)
    result_std_hate=hcat(time_list, hate_ratio_std)

    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_hate)
    writedlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_hate)
    
end

function g_ratio_evolution_NWOM_Sc_n_(type,time,num_of_individuals,beta,gamma,delta,alpha,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])
    
    #total_links=length(list_of_neighbors)/2
    
    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors
    
    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
        
    end

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    #dictionary_of_density_links=Dict{Float64, Vector{Float64}}()
    dictionary_of_haters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time
        dictionary_of_intervals[i]=[]
        #dictionary_of_density_links[i]=[]
        dictionary_of_new_adopters[i]=[]
        dictionary_of_haters[i]=[]
        push!(t_inter_arr, i)
    end


    ################################## MAIN LOOP ################################

    
    for _ in 1:cycle
        posibilities_array=[i for i in 1:2*num_of_individuals] # possible transitions for the non adopters (2 times the number of non adopters)

        t=0.0
        penetration_ratio=0.0
        hate_ratio=0.0
        status_list=zeros(num_of_individuals)

        w_array_adopt=ones(num_of_individuals).*gamma*(1-delta) # positive adoption
        w_array_hates=ones(num_of_individuals).*gamma*delta     # negative adoption
        w_arrays=vcat(w_array_adopt,w_array_hates)
        W=num_of_individuals*gamma
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0

        #density_links=0.0
        while t<time && !isempty(posibilities_array)

            t_n=-log(rand())/W
            t +=t_n
            r=rand()*W
            a=0
            adopter=posibilities_array[1]

            for n in posibilities_array

                a+=w_arrays[n]

                if r<a
                    adopter=n
                    break
                end
            end

            if adopter>num_of_individuals # the adopter hates the technology
                hate_ratio+=1/num_of_individuals
                adopter_2=Int(adopter-num_of_individuals) # position of the possibility of happy adoption 
                status_list[adopter_2]=-1.0

                neighbors=dictionary_of_neighbors[adopter_2] # neighbors of the adopter 
                status_of_neighborhood=[status_list[i] for i in dictionary_of_neighbors[adopter_2]] # state of the neighbors of the adopter

                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of unhappy adoption

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0

            else # the adopter is happy with the technology
                status_list[adopter]=1.0

                neighbors=dictionary_of_neighbors[adopter] # neighbors of the adopter
                status_of_neighborhood=[status_list[i] for i in dictionary_of_neighbors[adopter]] # state of the neighbors of the adopter

                adopter_2=Int(adopter+num_of_individuals) # position of the possibility of unhappy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of unhappy adoption

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0
            end

            penetration_ratio+=1/num_of_individuals

            counter=0
            for neighbor in neighbors

                neighbor_2=Int(neighbor+num_of_individuals) # position of the neighbor if it hates 
                counter+=1
                status_neighbor=status_of_neighborhood[counter]
                if status_neighbor==0.0 # if the neighbor is not infected we update the probability 

                    status_of_second_neighbors=[status_list[i] for i in dictionary_of_neighbors[neighbor]]                # length=k
                    positive_density=count(x -> x == 1.0, status_of_second_neighbors)/length(status_of_second_neighbors)  # n_+/k
                    negative_density=count(x -> x == -1.0, status_of_second_neighbors)/length(status_of_second_neighbors) # n_-/k

                    w_arrays[neighbor]=  (gamma+beta*positive_density)*(1-negative_density)^alpha*(1-delta)  # prob of positive adoption 
                    w_arrays[neighbor_2]=(gamma+beta*positive_density)*(1-negative_density)^alpha*delta      # prob of negative adoption 
                end
            end
            
            W=sum(w_arrays)
            
            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                push!(dictionary_of_haters[t_inter], hate_ratio)
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
    hate_mean=[]
    hate_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        #push!(density_links_mean, mean(dictionary_of_density_links[i]))
        #push!(density_links_std, std(dictionary_of_density_links[i]))
        push!(hate_mean, mean(dictionary_of_haters[i]))
        push!(hate_std, std(dictionary_of_haters[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)

    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_3)

    #result_density_links=hcat(t_inter_arr, density_links_mean)
    #writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_density_links)
    
    #result_std_density=hcat(t_inter_arr, density_links_std)
    #writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_density)

    result_hate=hcat(t_inter_arr, hate_mean)
    result_std_hate=hcat(t_inter_arr, hate_std)

    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_hate)
    writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_hate)
end

function f_ratio_evolution_NWOM_Sc_n_(type,time,num_of_individuals,beta,gamma,delta,alpha,cycle,degree_val)
    
    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])

    #total_links=length(list_of_neighbors)/2

    dictionary_of_neighbors=Dict{Int, Vector{Int}}() # empty dictionary, the keys are the individuals and the values are the list of neighbors

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighbors[i]=neighbors
    end

    ################################## PARAMETERS ################################

    K=100

    dictionary_of_intervals=Dict{Float64, Vector{Float64}}()
    dictionary_of_new_adopters=Dict{Float64, Vector{Float64}}()
    #dictionary_of_density_links=Dict{Float64, Vector{Float64}}()
    dictionary_of_haters=Dict{Float64, Vector{Float64}}()
    t_inter_arr=[]

    for i in 0:time/K:time
        dictionary_of_intervals[i]=[]
        #dictionary_of_density_links[i]=[]
        dictionary_of_new_adopters[i]=[]
        dictionary_of_haters[i]=[]
        push!(t_inter_arr, i)
    end

    ################################## MAIN LOOP ################################
    
    for _ in 1:cycle

        posibilities_array=[i for i in 1:2*num_of_individuals] # possible transitions for the non adopters (2 times the number of non adopters)

        t=0.0
        penetration_ratio=0.0
        hate_ratio=0.0
        status_list=zeros(num_of_individuals)

        w_array_adopt=ones(num_of_individuals).*gamma*(1-delta) # positive adoption
        w_array_hates=ones(num_of_individuals).*gamma*delta     # negative adoption
        w_arrays=vcat(w_array_adopt,w_array_hates)
        k_inter=1
        t_inter=t_inter_arr[k_inter]
        old_adopters=0.0


        while t<time && !isempty(posibilities_array)

            t_n=10^8
            adopter=0

            for n in posibilities_array

                local_rate=w_arrays[n]
                t_u=-log(rand())/local_rate
                if t_u<t_n
                    t_n=t_u
                    adopter=n
                end
            end

            if adopter>num_of_individuals # the adopter hates the technology

                hate_ratio+=1/num_of_individuals
                adopter_2=Int(adopter-num_of_individuals) # position of the possibility of happy adoption 

                neighbors=dictionary_of_neighbors[adopter_2] # neighbors of the adopter 
                status_of_neighborhood=[status_list[i] for i in dictionary_of_neighbors[adopter_2]] # state of the neighbors of the adopter

                status_list[adopter_2]=-1.0

                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of unhappy adoption

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0

            else # the adopter is happy with the technology
                status_list[adopter]=1.0

                neighbors=dictionary_of_neighbors[adopter] # neighbors of the adopter
                status_of_neighborhood=[status_list[i] for i in dictionary_of_neighbors[adopter]] # state of the neighbors of the adopter

                adopter_2=Int(adopter+num_of_individuals) # position of the possibility of unhappy adoption 
                posibilities_array=filter(!=(adopter),posibilities_array) # we erase the posibility of happy adoption 
                posibilities_array=filter(!=(adopter_2),posibilities_array) # we erase the posibility of unhappy adoption

                w_arrays[adopter]=0.0 # we set the probability of adoption to 0
                w_arrays[adopter_2]=0.0
            end

            penetration_ratio += 1/num_of_individuals
            t +=t_n

            counter=0
            for neighbor in neighbors

                neighbor_2=Int(neighbor+num_of_individuals)
                counter+=1
                status_neighbor=status_of_neighborhood[counter]
                if status_neighbor==0.0 # if the neighbor is not infected we update the probability 

                    status_of_second_neighbors=[status_list[i] for i in dictionary_of_neighbors[neighbor]]                # length=k
                    positive_density=count(x -> x == 1.0, status_of_second_neighbors)/length(status_of_second_neighbors)  # n_+/k
                    negative_density=count(x -> x == -1.0, status_of_second_neighbors)/length(status_of_second_neighbors) # n_-/k

                    w_arrays[neighbor]=  (gamma+beta*positive_density)*(1-negative_density)^alpha*(1-delta)  # prob of positive adoption 
                    w_arrays[neighbor_2]=(gamma+beta*positive_density)*(1-negative_density)^alpha*delta      # prob of negative adoption 
                end
            end
            
            if t>=t_inter
                push!(dictionary_of_intervals[t_inter], penetration_ratio)
                push!(dictionary_of_new_adopters[t_inter], penetration_ratio-old_adopters)
                push!(dictionary_of_haters[t_inter], hate_ratio)
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
    hate_mean=[]
    hate_std=[]
    for i in t_inter_arr
        push!(penetration_mean, mean(dictionary_of_intervals[i]))
        push!(penetration_std, std(dictionary_of_intervals[i]))
        push!(new_adopters_mean, mean(dictionary_of_new_adopters[i]))
        #push!(density_links_mean, mean(dictionary_of_density_links[i]))
        #push!(density_links_std, std(dictionary_of_density_links[i]))
        push!(hate_mean, mean(dictionary_of_haters[i]))
        push!(hate_std, std(dictionary_of_haters[i]))
    end

    results=hcat(t_inter_arr, penetration_mean)

    writedlm("NWOMBass_data//f_ratio_evolution_NWOM_Sc_n_-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", results)

    result_2=hcat(t_inter_arr, new_adopters_mean.*(K/time))
    writedlm("NWOMBass_data//f_ratio_evolution_NWOM_Sc_n_-new_adopt-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_2)

    result_3=hcat(t_inter_arr, penetration_std)
    writedlm("NWOMBass_data//f_ratio_evolution_NWOM_Sc_n_-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_3)

    #result_density_links=hcat(t_inter_arr, density_links_mean)
    #writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_density_links)
    
    #result_std_density=hcat(t_inter_arr, density_links_std)
    #writedlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std_density_links-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_density)

    result_hate=hcat(t_inter_arr, hate_mean)
    result_std_hate=hcat(t_inter_arr, hate_std)

    writedlm("NWOMBass_data//f_ratio_evolution_NWOM_Sc_n_-hate-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_hate)
    writedlm("NWOMBass_data//f_ratio_evolution_NWOM_Sc_n_-hate-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt", result_std_hate)

end

##############################

gamma=0.0004 #0.0004 
beta=1.088 #1.088 
alpha=2 #0.01 # alpha >= 1
delta=0.3
cycles=10^3
Num=10^4
time_g=30 
time=time_g*Num*(gamma+beta)

#= TYPES OF NETWORKS:

-square || 1
-erdos  ||0.1
-scale_free ||5
-z  ||5; 4

=#

#=
println("Discrete time algorithm")
@btime s_ratio_evolution_NWOM_Sc_n_("square",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 1)
println("First reaction algorithm")
@btime f_ratio_evolution_NWOM_Sc_n_("square",time_g,Num,beta,gamma,delta,alpha,cycles, 1)
println("Residence time algorithm")
@btime g_ratio_evolution_NWOM_Sc_n_("square",time_g,Num,beta,gamma,delta,alpha,cycles, 1)
# =#

#######################################################################

function plotting_adopters_MF_All_to_all(N,beta,gamma,delta,alpha,time)

    #result_1=readdlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-All_to_all-$N-$beta-$gamma-$delta-$alpha.txt")
    result_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-All_to_all-$N-$beta-$gamma-$delta-$alpha.txt")
    std_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std-All_to_all-$N-$beta-$gamma-$delta-$alpha.txt")

    mf_result=readdlm("NWOMBass_data//MF_NWOM_adopters_$beta-$gamma-$delta-$alpha.txt")
    
    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2,xlabel="Time",ylabel="Density of total adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//NWOM//Adopters//Adopters_all_to_all-$N-$beta-$gamma-$delta-$alpha.png")

end

function plotting_haters_MF_All_to_all(N,beta,gamma,delta,alpha,time)


    #result_1=readdlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-All_to_all-$N-$beta-$gamma-$delta-$alpha.txt")
    result_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-All_to_all-$N-$beta-$gamma-$delta-$alpha.txt")
    std_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-std-All_to_all-$Num-$beta-$gamma-$delta-$alpha.txt")

    mf_result=readdlm("NWOMBass_data//MF_hate_$beta-$gamma-$delta-$alpha.txt")

    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2, grid=false,xlabel="Time",ylabel="Density of negative adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//NWOM//Haters//Adopters_all_to_all-hate-$N-$beta-$gamma-$delta-$alpha.png")
end

#######################################################################

function plotting_adopters(type,N,beta,gamma,delta,alpha,degree_val,time)


    #result_1=readdlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")
    result_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")
    #std_1=readdlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-std-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")
    std_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-std-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")

    mf_result=readdlm("NWOMBass_data//MF_NWOM_adopters_$beta-$gamma-$delta-$alpha.txt")
        
    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2, grid=false,xlabel="Time",ylabel="Density of total adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//NWOM//Adopters//Adopters_$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.png")
end

function plotting_haters(type,N,beta,gamma,delta,alpha, degree_val,time)

    #result_1=readdlm("NWOMBass_data//s_ratio_evolution_NWOM_Sc_n_-hate-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")
    result_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-$type-$N-$beta-$gamma-$delta-$alpha-$degree_val.txt")
    std_2=readdlm("NWOMBass_data//g_ratio_evolution_NWOM_Sc_n_-hate-std-$type-$Num-$beta-$gamma-$delta-$alpha-$degree_val.txt")

    mf_result=readdlm("NWOMBass_data//MF_hate_$beta-$gamma-$delta-$alpha.txt")

    x_data=result_2[:,1]
    y_data=result_2[:,2]
    y_std=std_2[:,2]
    n = div(length(x_data), 10)
    x_subset = x_data[1:n:end]
    y_subset = y_data[1:n:end]
    y_std_subset = y_std[1:n:end]

    plot(result_2[:,1],result_2[:,2], label="Residence time algorithm", linewidth=4, linecolor=:magenta2, grid=false,xlabel="Time",ylabel="Density of negative adopters")
    plot!(mf_result[:,1],mf_result[:,2], label="MF approximation", linestyle=:dot, linewidth=4, linecolor=:green)
    plot!(x_subset, y_subset, yerr=y_std_subset, seriestype=:scatter, label="", markercolor=:magenta2)
    xlims!(0,time)
    savefig("Figures//NWOM//Haters//Adopters_$type-hate-$N-$beta-$gamma-$delta-$alpha-$degree_val.png")
end

#######################################################################

#=
s_ratio_evolution_NWOM_Sc_n_All_to_all(Int(round(time)),Num,beta,gamma,delta,alpha,cycles)
g_ratio_evolution_NWOM_Sc_n_All_to_all(time_g,Num,beta,gamma,delta,alpha,cycles)

s_ratio_evolution_NWOM_Sc_n_("square",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 1)
g_ratio_evolution_NWOM_Sc_n_("square",time_g,Num,beta,gamma,delta,alpha,cycles, 1)

s_ratio_evolution_NWOM_Sc_n_("scale_free",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 5)
g_ratio_evolution_NWOM_Sc_n_("scale_free",time_g,Num,beta,gamma,delta,alpha,cycles, 5)

s_ratio_evolution_NWOM_Sc_n_("z",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 5)
g_ratio_evolution_NWOM_Sc_n_("z",time_g,Num,beta,gamma,delta,alpha,cycles, 5)

s_ratio_evolution_NWOM_Sc_n_("z",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 4)
g_ratio_evolution_NWOM_Sc_n_("z",time_g,Num,beta,gamma,delta,alpha,cycles, 4)

s_ratio_evolution_NWOM_Sc_n_("erdos",Int(round(time)),Num,beta,gamma,delta,alpha,cycles, 0.1)
g_ratio_evolution_NWOM_Sc_n_("erdos",time_g,Num,beta,gamma,delta,alpha,cycles, 0.1)
# =#

# =
plotting_haters_MF_All_to_all(Num,beta,gamma,delta,alpha,time_g)
plotting_adopters_MF_All_to_all(Num,beta,gamma,delta,alpha,time_g)

plotting_adopters("square",Num,beta,gamma,delta,alpha,1,time_g)
plotting_haters("square",Num,beta,gamma,delta,alpha,1,time_g)

plotting_adopters("scale_free",Num,beta,gamma,delta,alpha,5,time_g)
plotting_haters("scale_free",Num,beta,gamma,delta,alpha,5,time_g)

plotting_adopters("z",Num,beta,gamma,delta,alpha,5,time_g)
plotting_haters("z",Num,beta,gamma,delta,alpha,5,time_g)

plotting_adopters("z",Num,beta,gamma,delta,alpha,4,time_g)
plotting_haters("z",Num,beta,gamma,delta,alpha,4,time_g)

plotting_adopters("erdos",Num,beta,gamma,delta,alpha,0.1,time_g)
plotting_haters("erdos",Num,beta,gamma,delta,alpha,0.1,time_g)
# =#