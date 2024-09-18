using DifferentialEquations, Plots, LaTeXStrings, Statistics, FindPeaks1D, DelimitedFiles

function unzip(tuples)
    a, b = [], []
    for (x, y) in tuples
        push!(a, x)
        push!(b, y)
    end
    return a, b
end

function count_decimals(x::Float64)
    # Convert the float to a string
    str = string(x)
    
    # Find the position of the decimal point
    dot_position = findfirst(==('.'), str)
    
    # If there is no decimal point, return 0
    if dot_position === nothing
        return 0
    end
    
    # Count the number of digits after the decimal point
    return length(str) - dot_position
end

function time_respect_parameters(beta, gamma,rho)
    t=1/(beta+gamma)*log((1+beta*rho/gamma)/(1-rho))
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end

function k_min_k_max(type,num_of_individuals,degree_val)

    pointer_to_neighbors=Int.(readdlm("networks_data//$type-pointer_$num_of_individuals-$degree_val.txt")[:,1])
    list_of_neighbors=Int.(readdlm("networks_data//$type-list_$num_of_individuals-$degree_val.txt")[:,1])

    dictionary_of_neighborhood=Dict{Int64, Vector{Int64}}()

    for i in 1:num_of_individuals
        neighbors=list_of_neighbors[pointer_to_neighbors[i]:pointer_to_neighbors[i+1]-1]
        dictionary_of_neighborhood[i]=neighbors
    end

    dictionary_of_degrees=Dict{Int64, Int64}()

    k_min=10^4
    k_max=0
    for i in 1:num_of_individuals
        k=length(dictionary_of_neighborhood[i])
        if k<k_min
            k_min=k
        end
        if k>k_max
            k_max=k
        end
        if k in keys(dictionary_of_degrees)
            dictionary_of_degrees[k]+=1
        else
            dictionary_of_degrees[k]=1
        end
    end

    writedlm("ratio_data/k_min_max-$type-$num_of_individuals-$degree_val.txt", [k_min, k_max])
    writedlm("ratio_data/k_degrees-$type-$num_of_individuals-$degree_val.txt", collect(values(dictionary_of_degrees)))
    writedlm("ratio_data/k_degrees_keys-$type-$num_of_individuals-$degree_val.txt", collect(keys(dictionary_of_degrees)))
end

function beta_s(s_array,gamma,beta,p_k,degree_array)
    numerator=0
    denominator=0
    j=0
    for k in degree_array
        j+=1
        for m in 1:k+1
            numerator+=p_k[j]*(k-m+1)*(gamma+beta*(m-1)/k)*s_array[k][m] # s_array is a dict 
            denominator+=p_k[j]*(k-m+1)*s_array[k][m]
        end
    end
    

    if denominator==0
        return 0
    else
        return numerator/denominator
    end
end

function AME_equations(type, beta,gamma,num_of_individuals, degree_val,time,dt) # Euler's method

    number_of_adopters_with_degree_k=readdlm("networks_data/k_degrees-$type-$num_of_individuals-$degree_val.txt")[:,1]
    degree_array=readdlm("networks_data/k_degrees_keys-$type-$num_of_individuals-$degree_val.txt")[:,1]

    number_of_adopters_with_degree_k=Int.(number_of_adopters_with_degree_k)
    degree_array=Int.(degree_array)

    # dictionary with the state of individuals given K over time and initial values

    fraction_non_adopters_k=Dict{Int64, Vector{Float64}}()

    for k in degree_array
        fraction_non_adopters_k[k]=zeros(k+1)  
    end


    #= fraction of k-degree nodes that is susceptible at time t and have m
    infected neighbors. 
    
    s_{k,0}=k/N; s_{k,m!=0}=0
    =#

    # to obtain the distribution p_k we need to zip together degree_array and number_of_adopters_with_degree_k

    combined = collect(zip(degree_array, number_of_adopters_with_degree_k))

    sorted_combined = sort(combined, by = x -> x[1])

    # and the unzip them

    degree_array, number_of_adopters_with_degree_k = unzip(sorted_combined)

    p_k=[]
    sum_degrees=0
    for k in degree_array
        position=findfirst(x->x==k, degree_array)
        fraction_non_adopters_k[k][1]=1 # fraction of non adopters of degree k that have no infected neighbors
        push!(p_k, number_of_adopters_with_degree_k[position]/num_of_individuals)
        sum_degrees+=k*number_of_adopters_with_degree_k[position]
    end

    #= dictionary that maps a times key to another dictionary. This second dictionary is the density of adopters
    with degree k and m=0,1,...,k in state 0. =#
    evolution_of_s_km=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()

    evolution_of_s_km[0.0]=[fraction_non_adopters_k]

    #=
    for t in dt:dt:time
        evolution_of_s_km[t]=[]
    end
    =#

    #=  We need a tripple nested loop
    -First loop(1/2): Loop over times
        In this loop we'll create a dictionary, to store the arrays with the state of the individuals for every k
        dictionary_states={k_min=>[s_{k_min,0}(t),s_{k_min,1}(t),...,s_{k_min,k}(t)],...,k_max=>[s_{k_max,0}(t),s_{k_max,1}(t),...,s_{k_max,k}(t)]}

    -Second loop(1/2): Loop over k
        In this loop we'll obtain the array of density of individuals in state 0 with degree k.
        We'll create two arrays, one that should not be updated in the next loop and another that should be.

    -Third loop(1/1): Loop over m 
        In this loop we'll obtain the denisity of individuals in state 0 at time t+dt for every value of m.
        We'll also update the array that should be updated with the new value of s_k_m(t+dt). 

    -Second loop (2/2): Loop over k
        Now, we'll update the array that should not be updated until the end.
        And we'll also create a new definition in the dictionary with the values of s_k_m(t+dt) for that value of k.
    
    -Third loop (2/2): Loop over times 
        Finally, we'll update evolution_of_s_km with the new values of s_k_m(t+dt) for every value of k.

    =#

    a_t=[0.0]              # initial values
    density_links_t=[0.0]

    for t in 0:dt:time-dt

        j=0
        summatory_pk=0
        summatory_density_links=0

        dictionary_states=Dict{Int64, Vector{Float64}}() # dictionary to store the state_array for every k
        beta_s_val=beta_s(evolution_of_s_km[t][1],gamma,beta,p_k,degree_array) #  β^s(t)

        for k in degree_array

            s_km_array=copy(evolution_of_s_km[t][1][k]) # s_k(t): fraction of nodes in state 0 with degree k (will be updated as time passes)
            s_updates=copy(evolution_of_s_km[t][1][k])
            sum_s_km=0 # array with  m*s_{k,m}:  0*s_{k,0};1*s_{k,1};2*s_{k,2};...;k*s_{k,k}
            
            for m in 1:k+1 # equivalent to sum_{m=0}^k

                s_km=s_km_array[m] # state of s_{k,m}(t)


                if m==1 #  state of s_{k,m-1}(t)
                    s_km_1=0
                else
                    s_km_1=s_km_array[m-1]
                end
                # s_{k,m}(t+dt) = s_{k,m} - F_{k,m} * s_{k,m}(t) * dt - β^s(t)*dt*(k-m)*s_{k,m}(t)+β^s*dt*(k-m+1)*s_{k,m-1}(t)
                s_km_dt=s_km - (gamma+beta*(m-1)/k) * s_km * dt -  beta_s_val * dt*(k-m+1)*s_km + beta_s_val * dt * (k-m+2)*s_km_1

                s_updates[m]=s_km_dt # update only the state of s_{k,m}(t+dt)
                sum_s_km+=(m-1)*s_km_dt
            end

    
            # all s_{k,m=0,1,...}(t+dt) have been updated, therefore we can update the complete array
            s_km_array=copy(s_updates) 

            dictionary_states[k]=s_km_array
            j+=1
            summatory_pk+=p_k[j]*sum(s_km_array)
            summatory_density_links+=p_k[j]*sum_s_km*num_of_individuals*2/sum_degrees
            
        end

        if t!=time
            evolution_of_s_km[round(t+dt, digits=number_digits)]=[dictionary_states]
        end

        #summatory_pk= \sum_k P_k \sum_{m=0}^k s_{k,m}(t_x)
        push!(a_t, 1-summatory_pk)
        push!(density_links_t,summatory_density_links)

        # we completely fill evolution_of_s_km[t+dt][1][k]. Now we must continue to the next t

    end



    #= Once we have complete the calculations for all k, we can obtain a(t).
        
        a(t)=1-\sum_k P_k \sum_m=0^k s_{k,m}(t)
        a(t_1)=1-
            [P_k_min*(s_{k_min,0}(t_1)+s_{k_min,1}(t_1)+...+s_{k,min,k}(t_1))
            +...+
            P_k_max*(s_{k_max,0}(t_1)+s_{k_max,1}(t_1)+...+s_{k_max,k}(t_1))]

        a_t=[a(0),a(dt),...,a(t)]=[1-\sum_k P_k \sum_{m=0}^k s_{k,m}(0),1-\sum_k P_k \sum_{m=0}^k s_{k,m}(dt),...,1-\sum_k P_k \sum_{m=0}^k s_{k,m}(t)]
    =#

    #=  We finally have in evolution_of_s_km[t][1][k] the density of adopters for every time, k and m in that order.
        Now we need to obtain \sum_{m=0}^k s_{k,m}(t)
    =#

    
    t_array=[t for t in 0:dt:time]

    #=
    a_t=[]
    for t in 0:dt:time # loop over times
        j=0
        summatory_pk=0

        s_km_t=evolution_of_s_km[t][1] # dictionary 

        for k in degree_array # loop over degrees
            j+=1
            summatory_pk+=p_k[j]*sum(s_km_t[k])
        end
        #summatory_pk= \sum_k P_k \sum_{m=0}^k s_{k,m}(t_x)
        push!(a_t, 1-summatory_pk)
    end
    =#

    result=hcat(t_array, a_t)
    writedlm("ratio_data/AME-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", result)
    result_2=hcat(t_array, density_links_t)
    writedlm("ratio_data/AME-density_links-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", result_2)
end

#= TYPES OF NETWORKS:

|NETWORK     || degree_val || dt    |
|------------||------------||-------|
|-square     || 1          || 0.01  |
|-erdos      ||0.1         || 0.001 |
|-scale_free ||5           || 0.001 |
|-z          ||5; 4        || 0.01  |


=#

Num=10^4

beta=1.088
gamma=0.0004
dt=0.01
time= time_respect_parameters(beta, gamma, 0.999)
number_digits=count_decimals(dt)

# =
#AME_equations("square",beta,gamma,Num,1,time,dt)
#AME_equations("z",beta,gamma,Num,4,time,dt)
#AME_equations("z",beta,gamma,Num,5,time,dt)
#AME_equations("scale_free",beta,gamma,Num,5,time,dt)
#AME_equations("All_to_all",beta,gamma,Num,"-",time,dt)
#AME_equations("erdos",beta,gamma,Num,0.1,time,dt)
# =#



