using Plots, LaTeXStrings, Statistics,DelimitedFiles,Serialization

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

function delta(k_1,k_2) # kronecker delta
    if k_1==k_2
        return 1.0
    else
        return 0.0
    end
end

function denominator_kq(degree_array,phi_dict)
    sum_denominator=0

    for k in degree_array
        for q in 1:k+1
            sum_denominator+=(k-q+1)*phi_dict[k][q]
        end
    end

    return sum_denominator
end

function denominator_q(degree_array,phi_dict)
    sum_denominator=0

    for k in degree_array
        for q in 1:k+1
            sum_denominator+=(q-1)*phi_dict[k][q]
        end
    end

    return sum_denominator
end
################################################ PROBABILITIES ############################################

function p00(k1,q1,phi0,denominator)
    numerator=(k1-q1)*phi0
    
    if denominator==0
        return 0.0
    else
        return numerator/denominator
    end
end

function p01(k1,q1,phi1,denominator)
    numerator=q1*phi1 #(k1-q1)*phi1 

    if denominator==0
        return 0.0
    else
        return numerator/denominator
    end
end

################################################ RATES ################################################
function beta_s(phi_dict,gamma,beta,degree_array,denominator)
    numerator=0
    for k in degree_array
        for q in 1:k+1
            numerator+=(k-q+1)*(gamma+beta*(q-1)/k)*phi_dict[k][q]
        end
    end

    if denominator==0
        return 0.0
    else
        return numerator/denominator
    end
end

function beta_i(phi_dict,gamma,beta,degree_array,denominator)
    numerator=0
    for k in degree_array
        for q in 1:k+1
            numerator+=(q-1)*(gamma+beta*(q-1)/k)*phi_dict[k][q]
        end
    end
    
    if denominator==0
        return 0.0
    else
        return numerator/denominator
    end
    
end

function beta_nss(phi_dict,gamma,beta,degree_array)
    numerator=0
    for k in degree_array     
        for q in 1:k+1
            numerator+=(k-q+1)*(k-q)*(gamma+beta*(q-1)/k)*phi_dict[k][q]
        end
    end
    
    return numerator
end

function beta_ns(phi_dict,gamma,beta,degree_array)
    numerator=0
    
    for k in degree_array
        for q in 1:k+1
            numerator+=(k-q+1)*(gamma+beta*(q-1)/k)*phi_dict[k][q]
        end
    end
    
    return numerator
end

function beta_nsi(phi0_dict,gamma,beta,degree_array)
    numerator=0
    
    for k in degree_array
        for q in 1:k+1
            numerator+=(q-1)*(k-q+1)*(gamma+beta*(q-1)/k)*phi0_dict[k][q]
        end
    end
    
    return numerator
end

function beta_nii(phi0_dict,gamma,beta,degree_array)
    numerator=0
    
    for k in degree_array       
        for q in 1:k+1
            numerator+=(q-1)*(q-2)*(gamma+beta*(q-1)/k)*phi0_dict[k][q]
        end
    end
    
    return numerator
end

function beta_ni(phi0_dict,gamma,beta,degree_array)
    numerator=0
    
    for k in degree_array
        for q in 1:k+1
            numerator+=(q-1)*(gamma+beta*(q-1)/k)*phi0_dict[k][q]
        end
    end
    
    return numerator
end

################################################ DERIVATIVES ############################################

function first_beta_s(gamma,beta,k1,q1,beta_s_val, denominator_kq)

    numerator=(k1-q1)*((gamma+beta*q1/k1)-beta_s_val)

    if denominator_kq==0
        return 0.0
    else
        return numerator/denominator_kq
    end
end

function first_beta_i(gamma,beta,k1,q1,beta_i_val,denominator_q)

    numerator=q1*((gamma+beta*q1/k1)-beta_i_val)

    if denominator_q==0
        return 0.0
    else
        return numerator/denominator_q
    end
end

function second_beta_s(gamma,beta,k1,q1,k2,q2,beta_s_val,denominator)
    numerator=(k1-q1)*(k2-q2)*(2*beta_s_val-(gamma+beta*q1/k1)-(gamma+beta*q2/k2))

    if denominator==0
        return 0.0
    else
        return numerator/denominator^2
    end
end

function second_beta_i(gamma,beta,k1,q1,k2,q2,beta_i_val,denominator)
    numerator=q1*q2*(2*beta_i_val-(gamma+beta*q1/k1)-(gamma+beta*q2/k2))

    if denominator==0
        return 0.0
    else
        return numerator/denominator^2
    end
end

################################################ MATRIX B ############################################

function b00(k1,k2,q1,q2,phi0,phi0_1,beta_s_val,first_beta_s_val,gamma,beta)
    return(
        delta(k1,k2)*(delta(q1,q2)*(gamma+beta*q1/k1+(k1-q1)*beta_s_val)-
        delta(q1-1,q2)*(k1-q1+1)*beta_s_val)+
        ((k1-q1)*phi0-(k1-q1+1)*phi0_1)*first_beta_s_val
        )
end

function b10(k1,k2,q1,q2,phi1,phi1_1,first_beta_i_val,gamma,beta)
    return(
        -delta(k1,k2)*delta(q1,q2)*(gamma+beta*q1/k1)+((k1-q1)*phi1-(k1-q1+1)*phi1_1)*first_beta_i_val)
end

function  b11(k1,k2,q1,q2,beta_i_val)
    return(
        delta(k1,k2)*(delta(q1,q2)*(k1-q1)-
        delta(q1-1,q2)*(k1-q1+1))*beta_i_val)
end

################################################ MATRIX G ############################################

function g00(k1,k2,q1,q2,phi01,phi02,p00_dict,beta_nss_val,beta_ns_val,gamma,beta)

    p00_k1_q1=p00_dict[k1][q1+1]
    p00_k2_q2=p00_dict[k2][q2+1]

    if q1==0
        p00_k1_q1_1=0
    else
        p00_k1_q1_1=p00_dict[k1][q1]
    end

    if q2==0
        p00_k2_q2_1=0
    else
        p00_k2_q2_1=p00_dict[k2][q2]
    end

    return(
        phi01*(gamma+beta*q1/k1)*delta(k1,k2)*delta(q1,q2)-
        phi01*(gamma+beta*q1/k1)*(k1-q1)*(-p00_k2_q2+p00_k2_q2_1)-
        phi02*(gamma+beta*q2/k2)*(k2-q2)*(-p00_k1_q1+p00_k1_q1_1)+
        beta_nss_val*(p00_k1_q1*p00_k2_q2-
        p00_k1_q1_1*p00_k2_q2-p00_k1_q1*p00_k2_q2_1+
        p00_k1_q1_1*p00_k2_q2_1)+
        beta_ns_val*delta(k1,k2)*(delta(q1,q2)*p00_k1_q1-
        delta(q1-1,q2)*p00_k1_q1_1-delta(q1+1,q2)*p00_k1_q1+
        delta(q1,q2)*p00_k1_q1_1)
    )
end

function g01(k1,k2,q1,q2,phi01,phi02,p00_dict,p01_dict,beta_nsi_val,gamma,beta)

    p00_k1_q1=p00_dict[k1][q1+1]
    p01_k2_q2=p01_dict[k2][q2+1]

    if q1==0
        p00_k1_q1_1=0
    else
        p00_k1_q1_1=p00_dict[k1][q1]
    end

    if q2==0
        p01_k2_q2_1=0
    else
        p01_k2_q2_1=p01_dict[k2][q2]
    end
    
    return(
        -phi01*(gamma+beta*q1/k1)*delta(k1,k2)*delta(q1,q2)-
        phi01*(gamma+beta*q1/k1)*q1*(-p01_k2_q2+p01_k2_q2_1)+
        phi02*(gamma+beta*q2/k2)*(k2-q2)*(-p00_k1_q1+p00_k1_q1_1)+
        beta_nsi_val*(p00_k1_q1*p01_k2_q2-
        p00_k1_q1*p01_k2_q2_1-p00_k1_q1_1*p01_k2_q2+
        p00_k1_q1_1*p01_k2_q2_1)
    )
end

function g11(k1,k2,q1,q2,phi01,phi02,p01_dict,beta_nii_val,beta_ni_val,gamma,beta)

    p01_k1_q1=p01_dict[k1][q1+1]
    p01_k2_q2=p01_dict[k2][q2+1]

    if q1==0
        p01_k1_q1_1=0
    else
        p01_k1_q1_1=p01_dict[k1][q1]
    end

    if q2==0
        p01_k2_q2_1=0
    else
        p01_k2_q2_1=p01_dict[k2][q2]
    end

    return(
        phi01*(gamma+beta*q1/k1)*delta(k1,k2)*delta(q1,q2)+
        phi01*(gamma+beta*q1/k1)*q1*(-p01_k2_q2+p01_k2_q2_1)+
        phi02*(gamma+beta*q2/k2)*q2*(-p01_k1_q1+p01_k1_q1_1)+
        beta_nii_val*(p01_k1_q1*p01_k2_q2-
        p01_k1_q1_1*p01_k2_q2-p01_k1_q1*p01_k2_q2_1+p01_k1_q1_1*p01_k2_q2_1)+
        beta_ni_val*delta(k1,k2)*(delta(q1,q2)*p01_k1_q1-
        delta(q1-1,q2)*p01_k1_q1_1-delta(q1+1,q2)*p01_k1_q1+delta(q1,q2)*p01_k1_q1_1)
    )
    
end

################################################ MATRIX GAMMA ################################################

function Phi0_00(k1,k2,k3,q1,q2,q3,phi01,phi01_1,first_beta_s_val_k2q2,first_beta_s_val_k3q3,second_beta_s_val)

    return(
        delta(k1,k2)*(-delta(q1,q2)*(k1-q1)+delta(q1-1,q2)*(k1-q1+1))*first_beta_s_val_k3q3+
        delta(k1,k3)*(-delta(q1,q3)*(k1-q1)+delta(q1-1,q3)*(k1-q1+1))*first_beta_s_val_k2q2+
        (-(k1-q1)*phi01+(k1-q1+1)*phi01_1)*second_beta_s_val
    )
    
end

function Phi1_00(k1,q1,phi1,phi1_1,second_beta_i_val)

    return(
        (-(k1-q1)*phi1+(k1-q1+1)*phi1_1)*second_beta_i_val
    )
    
end

function Phi1_01(k1,k3,q1,q3,first_beta_i_val)

    return(
        delta(k1,k3)*(-delta(q1,q3)*(k1-q1)+delta(q1-1,q3)*(k1-q1+1))*first_beta_i_val
    )
    
end

#############################################################################################################

function SAME_phi_equations(type, beta,gamma,num_of_individuals, degree_val,time,dt) # Euler's method

    number_of_adopters_with_degree_k=readdlm("ratio_data/k_degrees-$type-$num_of_individuals-$degree_val.txt")[:,1]
    degree_array=readdlm("ratio_data/k_degrees_keys-$type-$num_of_individuals-$degree_val.txt")[:,1]

    number_of_adopters_with_degree_k=Int.(number_of_adopters_with_degree_k)
    degree_array=Int.(degree_array)

    # dictionary with the state of individuals given K over time and initial values

    fraction_non_adopters_k=Dict{Int64, Vector{Float64}}()
    fraction_adopters_k=Dict{Int64, Vector{Float64}}()


    for k in degree_array
        fraction_non_adopters_k[k]=zeros(k+1)  
        fraction_adopters_k[k]=zeros(k+1)  
    end


    #= fraction of k-degree nodes that is susceptible at time t and have m
    infected neighbors. 
    
    s_{k,0}=k/N; s_{k,m!=0}=0
    =#
    
    combined = collect(zip(degree_array, number_of_adopters_with_degree_k))

    sorted_combined = sort(combined, by = x -> x[1])

    degree_array, number_of_adopters_with_degree_k = unzip(sorted_combined)

    for k in degree_array
        position=findfirst(x->x==k, degree_array)
        # fraction of non adopters of degree k that have no infected neighbors:
        fraction_non_adopters_k[k][1]= number_of_adopters_with_degree_k[position]/num_of_individuals 
    end

    #= dictionary that maps a times key to another dictionary. This second dictionary is the density of adopters
    with degree k and m=0,1,...,k in state 0. =#
    evolution_of_phi0_kq=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    evolution_of_phi1_kq=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()

    evolution_of_phi0_kq[0.0]=[fraction_non_adopters_k] # initial condition \phi_{0}
    evolution_of_phi1_kq[0.0]=[fraction_adopters_k]     # initial condition \phi_{1}

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

    beta_s_array=[]
    beta_i_array=[]
    beta_nss_array=[]
    beta_ns_array=[]
    beta_nsi_array=[]
    beta_nii_array=[]
    beta_ni_array=[]

    denominator_kq_array=[] # \sum_k\sum_{q=0}^k (k-q)*phi_{0,k,q}
    denominator_kq_1_array=[] # \sum_k\sum_{q=0}^k (k-q)*phi_{1,k,q}
    denominator_q_array=[] # \sum_k\sum_{q=0}^k q*phi_{0,k,q}
    denominator_q_1_array=[] # \sum_k\sum_{q=0}^k q*phi_{1,k,q}

    
    for t in 0:dt:time
        t_dt=round(t+dt, digits=number_digits)
        dictionary_states_0=Dict{Int64, Vector{Float64}}() # dictionary to store the state_array for every phi0_k
        dictionary_states_1=Dict{Int64, Vector{Float64}}() # dictionary to store the state_array for every phi1_k
        
        denominator_kq_val=denominator_kq(degree_array, evolution_of_phi0_kq[t][1])
        denominator_kq_1_val=denominator_kq(degree_array, evolution_of_phi1_kq[t][1])
        denominator_q_val=denominator_q(degree_array, evolution_of_phi0_kq[t][1])
        denominator_q_1_val=denominator_q(degree_array, evolution_of_phi1_kq[t][1])
        beta_s_val=    beta_s(evolution_of_phi0_kq[t][1],gamma,beta,degree_array,denominator_kq_val) #  β^s(t)
        beta_i_val=    beta_i(evolution_of_phi0_kq[t][1],gamma,beta,degree_array,denominator_q_val) #  β^i(t)
        beta_nss_val=beta_nss(evolution_of_phi0_kq[t][1],gamma,beta,degree_array)
        beta_ns_val=  beta_ns(evolution_of_phi0_kq[t][1],gamma,beta,degree_array)
        beta_nsi_val=beta_nsi(evolution_of_phi0_kq[t][1],gamma,beta,degree_array)
        beta_nii_val=beta_nii(evolution_of_phi0_kq[t][1],gamma,beta,degree_array)
        beta_ni_val=  beta_ni(evolution_of_phi0_kq[t][1],gamma,beta,degree_array)

        push!(denominator_kq_array, denominator_kq_val)
        push!(denominator_kq_1_array, denominator_kq_1_val)
        push!(denominator_q_array, denominator_q_val)
        push!(denominator_q_1_array, denominator_q_1_val)
        push!(beta_s_array,beta_s_val)
        push!(beta_i_array,beta_i_val)
        push!(beta_nss_array,beta_nss_val)
        push!(beta_ns_array,beta_ns_val)
        push!(beta_nsi_array,beta_nsi_val)
        push!(beta_nii_array,beta_nii_val)
        push!(beta_ni_array,beta_ni_val)

        for k in degree_array

            #phi0_k(t): fraction of nodes in state 0 with degree k for all q's(will be updated as time passes)
            phi0_kq_array=copy(evolution_of_phi0_kq[t][1][k]) 
            phi0_updates=copy(evolution_of_phi0_kq[t][1][k])

            #phi1_k(t): fraction of nodes in state 0 with degree k for all q's(will be updated as time passes)
            phi1_kq_array=copy(evolution_of_phi1_kq[t][1][k]) 
            phi1_updates=copy(evolution_of_phi1_kq[t][1][k])


            for q in 1:k+1 # equivalent to sum_{q=0}^k

                phi0_kq=phi0_kq_array[q] # state of phi0_{k,q}(t)
                phi1_kq=phi1_kq_array[q] # state of phi1_{k,q}(t)

                if q==1 #  state of phi_{n,k,q-1}(t)
                    phi0_kq_1=0
                    phi1_kq_1=0
                else
                    phi0_kq_1=phi0_kq_array[q-1]
                    phi1_kq_1=phi1_kq_array[q-1]
                end

                phi0_kq_dt=phi0_kq - (gamma+beta*(q-1)/k) *phi0_kq * dt - beta_s_val * dt*(k-q+1)*phi0_kq + beta_s_val * dt * (k-q+2)*phi0_kq_1
                phi1_kq_dt=phi1_kq + (gamma+beta*(q-1)/k) *phi0_kq * dt - beta_i_val * dt*(k-q+1)*phi1_kq + beta_i_val * dt * (k-q+2)*phi1_kq_1

                phi0_updates[q]=phi0_kq_dt # update only the state of phi0_{k,m}
                phi1_updates[q]=phi1_kq_dt # update only the state of phi1_{k,m}

            end

    
            # all phi_{n,k,q=0,1,...} have been updated, therefore we can update the complete array
            phi0_kq_array=copy(phi0_updates)
            phi1_kq_array=copy(phi1_updates) 

            dictionary_states_0[k]=phi0_kq_array
            dictionary_states_1[k]=phi1_kq_array
            
        end

        evolution_of_phi0_kq[t_dt]=[dictionary_states_0]
        evolution_of_phi1_kq[t_dt]=[dictionary_states_1]

        # we completely fill phi_{n,k,q}[t+dt]. Now we must continue to the next t

    end

    open("ratio_data/phi0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, evolution_of_phi0_kq)
    end

    open("ratio_data/phi1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, evolution_of_phi1_kq)
    end

    writedlm("ratio_data/denominator_kq-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", denominator_kq_array)
    writedlm("ratio_data/denominator_kq_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", denominator_kq_1_array)
    writedlm("ratio_data/denominator_q-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", denominator_q_array)  
    writedlm("ratio_data/denominator_q_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", denominator_q_1_array)    
    writedlm("ratio_data/beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_s_array)
    writedlm("ratio_data/beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_i_array)
    writedlm("ratio_data/beta_nss-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_nss_array)
    writedlm("ratio_data/beta_ns-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_ns_array)
    writedlm("ratio_data/beta_nsi-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_nsi_array)
    writedlm("ratio_data/beta_nii-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_nii_array)
    writedlm("ratio_data/beta_ni-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", beta_ni_array)

    writedlm("ratio_data/degree_array-$type-$degree_val.txt", degree_array)
end

# function to compute the first derivative of β^s and β^i also with p_0(0,k,q) and p_0(1,k,q) over time
function SAME_first_bs_bi_and_p00_p01(type, beta,gamma,num_of_individuals, degree_val,time,dt)

    phi0_dictionary= open("ratio_data/phi0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    phi1_dictionary= open("ratio_data/phi1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    beta_s_array=readdlm("ratio_data/beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_i_array=readdlm("ratio_data/beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]

    denominator_kq_array=readdlm("ratio_data/denominator_kq-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    denominator_q_array=readdlm("ratio_data/denominator_q-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    denominator_kq_1_array=readdlm("ratio_data/denominator_kq_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    denominator_q_1_array=readdlm("ratio_data/denominator_q_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    

    # ∂b^s/∂ϕ_{0,k,q} ; ∂b^i/∂ϕ_{0,k,q} ; P00 ; P01 ; P10 ; P11  ###########################

    first_derivative_beta_s=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    first_derivative_beta_i=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    p00_dict_t=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    p01_dict_t=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()

    counter=1
    for t in 0:dt:time
        #t_dt=round(t, digits=3)
        dictionary_states_s=Dict{Int64, Vector{Float64}}() # dictionary to store first derivative b^s_{k,m}(t)
        dictionary_states_i=Dict{Int64, Vector{Float64}}() # dictionary to store first derivative b^i_{k,m}(t)

        dictionary_states_00=Dict{Int64, Vector{Float64}}() # dictionary to store p_0(0,k,q)(t)
        dictionary_states_01=Dict{Int64, Vector{Float64}}() # dictionary to store p_0(0,k,q)(t)
        
        beta_s_val=beta_s_array[counter] #  β^s(t)
        beta_i_val=beta_i_array[counter] #  β^i(t)
        denominator_kq_val=denominator_kq_array[counter]
        denominator_q_val=denominator_q_array[counter]
        denominator_q_1_val=denominator_q_1_array[counter]
        counter+=1

        for k in degree_array

            phi0_kq_array=copy(phi0_dictionary[t][1][k]) # array with ϕ_{0,k,_}(t) for all q
            phi1_kq_array=copy(phi1_dictionary[t][1][k]) # array with ϕ_{1,k,_}(t) for all q

            first_beta_s_array=[]
            first_beta_i_array=[]

            p00_array=[]
            p01_array=[]

            for q in 1:k+1

                phi0=phi0_kq_array[q] # ϕ_{0,k,q}
                phi1=phi1_kq_array[q] # ϕ_{1,k,q}
                
                push!(first_beta_s_array,first_beta_s(gamma,beta,k,q-1,beta_s_val, denominator_kq_val))
                push!(first_beta_i_array,first_beta_i(gamma,beta,k,q-1,beta_i_val, denominator_q_val))

                push!(p00_array, p00(k,q-1,phi0,denominator_kq_val))
                push!(p01_array, p01(k,q-1,phi1,denominator_q_1_val))
            end
            #= now we have the arrays with all the values of [∂b^s/∂ϕ_{0,k,q}](t) and [∂b^i/∂ϕ_{0,k,q}](t) for a given k and all q's
            we have to store those array in their respectives dictionaries =#
            dictionary_states_s[k]=first_beta_s_array
            dictionary_states_i[k]=first_beta_i_array

            # the same with P00, P01
            dictionary_states_00[k]=p00_array
            dictionary_states_01[k]=p01_array
        end

        #= once we have all the values of [∂b^s/∂ϕ_{0,k,q}](t) and [∂b^i/∂ϕ_{0,k,q}](t) for all k, 
        we store them inside their respectives dictionarties =#

        first_derivative_beta_s[t]=[dictionary_states_s]
        first_derivative_beta_i[t]=[dictionary_states_i]
        p00_dict_t[t]=[dictionary_states_00]
        p01_dict_t[t]=[dictionary_states_01]
    end

    open("ratio_data/p00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, p00_dict_t)
    end

    open("ratio_data/p01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, p01_dict_t)
    end

    open("ratio_data/first_derivative_beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, first_derivative_beta_s)
    end

    open("ratio_data/first_derivative_beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, first_derivative_beta_i)
    end

end

# function to compute the values of the B and G matrices over time
function SAME_B00_B10_B11_G00_G01_G10_G11(type, beta,gamma,num_of_individuals, degree_val,time,dt)

    phi0_dictionary= open("ratio_data/phi0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    phi1_dictionary= open("ratio_data/phi1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    beta_s_array=readdlm("ratio_data/beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_i_array=readdlm("ratio_data/beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_nss_array=readdlm("ratio_data/beta_nss-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_ns_array=readdlm("ratio_data/beta_ns-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_nsi_array=readdlm("ratio_data/beta_nsi-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_nii_array=readdlm("ratio_data/beta_nii-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_ni_array=readdlm("ratio_data/beta_ni-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]


    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    

    first_derivative_beta_s=open("ratio_data/first_derivative_beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    first_derivative_beta_i=open("ratio_data/first_derivative_beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    p00_dict_t=open("ratio_data/p00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    p01_dict_t=open("ratio_data/p01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    #                            t   , k , q , k', q'
    dict_Tuple_B00=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_B10=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_B11=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_G00=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_G01=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_G10=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_Tuple_G11=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()

    counter=1
    for t in 0:dt:time

        beta_s_val=beta_s_array[counter]     #  β^s(t)
        beta_i_val=beta_i_array[counter]     #  β^i(t)
        beta_nss_val=beta_nss_array[counter] # b^{nss}(t)
        beta_ns_val=beta_ns_array[counter]   # b^{ns}(t)
        beta_nsi_val=beta_nsi_array[counter] # b^{nsi}(t)
        beta_nii_val=beta_nii_array[counter] # b^{nii}(t)
        beta_ni_val=beta_ni_array[counter]   # b^{ni}(t)

        counter+=1

        p00_dict=p00_dict_t[t][1] # p_0(0,k,q)(t)
        p01_dict=p01_dict_t[t][1] # p_0(1,k,q)(t)

        for k in degree_array

            for q in 1:k+1

                
                phi0=copy(phi0_dictionary[t][1][k][q]) # ϕ_{0,k,q}
                phi1=copy(phi1_dictionary[t][1][k][q]) # ϕ_{1,k,q}
                if q==1
                    phi0_1=0
                    phi1_1=0
                else
                    phi0_1=copy(phi0_dictionary[t][1][k][q-1]) # ϕ_{0,k,q-1}
                    phi1_1=copy(phi1_dictionary[t][1][k][q-1]) # ϕ_{1,k,q-1}
                end

                for k_p in degree_array

                    for q_p in 1:k_p+1
                        first_beta_s_val=copy(first_derivative_beta_s[t][1][k_p][q_p]) # ∂b^s/∂ϕ_{0,k',q'}
                        first_beta_i_val=copy(first_derivative_beta_i[t][1][k_p][q_p]) # ∂b^i/∂ϕ_{0,k',q'}

                        phi0_p=copy(phi0_dictionary[t][1][k_p][q_p]) # ϕ_{0,k',q'}

                        val_B00=b00(k,k_p,q-1,q_p-1,phi0,phi0_1,beta_s_val,first_beta_s_val,gamma,beta)
                        val_B10=b10(k,k_p,q-1,q_p-1,phi1,phi1_1,first_beta_i_val,gamma,beta)
                        val_B11=b11(k,k_p,q-1,q_p-1,beta_i_val)

                        val_G00=g00(k,k_p,q-1,q_p-1,phi0,phi0_p,p00_dict,beta_nss_val,beta_ns_val,gamma,beta)
                        val_G11=g11(k,k_p,q-1,q_p-1,phi0,phi0_p,p01_dict,beta_nii_val,beta_ni_val,gamma,beta)
                        val_G01=g01(k,k_p,q-1,q_p-1,phi0,phi0_p,p00_dict,p01_dict,beta_nsi_val,gamma,beta)
                        val_G10=g01(k_p,k,q_p-1,q-1,phi0_p,phi0,p00_dict,p01_dict,beta_nsi_val,gamma,beta)
                        
                        tuple_parameters=(t,k,q,k_p,q_p)

                        dict_Tuple_B00[tuple_parameters]=val_B00
                        dict_Tuple_B10[tuple_parameters]=val_B10
                        dict_Tuple_B11[tuple_parameters]=val_B11

                        dict_Tuple_G00[tuple_parameters]=val_G00
                        dict_Tuple_G01[tuple_parameters]=val_G01
                        dict_Tuple_G10[tuple_parameters]=val_G10
                        dict_Tuple_G11[tuple_parameters]=val_G11

                    end
                end
            end
        end
    end
    # MATRIX B COMPONENTS ###
    open("ratio_data/B00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_B00)
    end
    open("ratio_data/B10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_B10)
    end
    open("ratio_data/B11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_B11)
    end
    # MATRIX G COMPONENTS ###
    open("ratio_data/G00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_G00)
    end
    open("ratio_data/G01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_G01)
    end
    open("ratio_data/G10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_G10)
    end
    open("ratio_data/G11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Tuple_G11)
    end
end

# function to compute the values of the C matrix over time
function SAME_C00_C01_C10_C11(type, beta,gamma,num_of_individuals, degree_val,time,dt)

    dictionary_B00= open("ratio_data/B00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_B10= open("ratio_data/B10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_B11= open("ratio_data/B11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_G00= open("ratio_data/G00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_G01= open("ratio_data/G01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_G10= open("ratio_data/G10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_G11= open("ratio_data/G11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    
    
    #                            t    , k , q , k', q'
    dict_tuple_C00=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_tuple_C01=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_tuple_C10=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    dict_tuple_C11=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()

    # initial conditions at t=0 C_{n,k,q;n',k',q'}(0)=0
    for k in degree_array
        for q in 1:k+1
            for k_p in degree_array
                for q_p in 1:k_p+1
                    tuple_param=(0,k,q,k_p,q_p)
                    dict_tuple_C00[tuple_param]=0.0
                    dict_tuple_C01[tuple_param]=0.0
                    dict_tuple_C10[tuple_param]=0.0
                    dict_tuple_C11[tuple_param]=0.0
                end
            end
        end
    end

    
    for t in 0:dt:time-dt
        t_dt=round(t+dt, digits=number_digits)
        for k in degree_array
            for q in 1:k+1
                for k_p in degree_array
                    for q_p in 1:k_p+1
                        tuple_params=(t,k,q,k_p,q_p) # (k, q; k', q')
                        
                        # =
                        sum_c00=0
                        sum_c01=0
                        sum_c10=0
                        sum_c11=0
                        for k_pp in degree_array
                            
                            for q_pp in 1:k_pp+1
                                tuple_params_ik=(t,k,q,k_pp,q_pp)     # (k, q ; k'', q'')
                                tuple_params_jk=(t,k_p,q_p,k_pp,q_pp) # (k',q'; k'', q'')

                                c00_ik=copy(dict_tuple_C00[tuple_params_ik])
                                c00_jk=copy(dict_tuple_C00[tuple_params_jk])
                                c01_ik=copy(dict_tuple_C01[tuple_params_ik])
                                c01_jk=copy(dict_tuple_C01[tuple_params_jk])
                                c10_ik=copy(dict_tuple_C10[tuple_params_ik])
                                c10_jk=copy(dict_tuple_C10[tuple_params_jk])
                                c11_ik=copy(dict_tuple_C11[tuple_params_ik])
                                c11_jk=copy(dict_tuple_C11[tuple_params_jk])

                                b00_ik=copy(dictionary_B00[tuple_params_ik])
                                b00_jk=copy(dictionary_B00[tuple_params_jk])
                                b10_ik=copy(dictionary_B10[tuple_params_ik])
                                b10_jk=copy(dictionary_B10[tuple_params_jk])
                                b11_ik=copy(dictionary_B11[tuple_params_ik])
                                b11_jk=copy(dictionary_B11[tuple_params_jk])

                                sum_c00+=(c00_jk*b00_ik+c00_ik*b00_jk)
                                sum_c01+=(c10_jk*b00_ik+c00_ik*b10_jk+c01_ik*b11_jk)
                                sum_c10+=(c00_jk*b10_ik+c10_ik*b00_jk+c01_jk*b11_ik)
                                sum_c11+=(c10_jk*b10_ik+c10_ik*b10_jk+c11_jk*b11_ik+c11_ik*b11_jk)

                            end
                        end
                        # =#
                        
                        c00_val=copy(dict_tuple_C00[tuple_params])
                        c01_val=copy(dict_tuple_C01[tuple_params])
                        c10_val=copy(dict_tuple_C10[tuple_params])
                        c11_val=copy(dict_tuple_C11[tuple_params])

                        #=
                        b00_val=copy(dictionary_B00[tuple_params])
                        b10_val=copy(dictionary_B10[tuple_params])
                        b11_val=copy(dictionary_B11[tuple_params])
                        # =#

                        #CB^T+BC ####
                        
                        #=
                        sum_c00=2*c00_val*b00_val
                        sum_c01=c00_val*b10_val+c01_val*b11_val+b00_val*c01_val
                        sum_c10=c10_val*b00_val+b10_val*c00_val+b11_val*c10_val
                        sum_c11=c10_val*b10_val+2*c11_val*b11_val+b10_val*c01_val
                        # =#
                        #CB+BC ####
                        #=
                        sum_c00=2*c00_val*b00_val+c01_val*b10_val
                        sum_c01=c01_val*b11_val+b00_val*c01_val
                        sum_c10=c10_val*b00_val+c11_val*b10_val+b10_val*c00_val+b11_val*c10_val
                        sum_c11=2*c11_val*b11_val+b10_val*c01_val
                        # =#
                        #CB^T+BC^T ####
                        #=
                        sum_c00=2*c00_val*b00_val
                        sum_c01=c00_val*b10_val+c01_val*b11_val+b00_val*c10_val
                        sum_c10=c10_val*b00_val+b10_val*c00_val+b11_val*c01_val
                        sum_c11=2*c10_val*b10_val+2*c11_val*b11_val
                        # =#

                        g00_val=copy(dictionary_G00[tuple_params]) # G_{0,k,q;0,k',q'}(t)
                        g01_val=copy(dictionary_G01[tuple_params]) # G_{0,k,q;1,k',q'}(t)
                        g10_val=copy(dictionary_G10[tuple_params]) # G_{0,k,q;1,k',q'}(t)
                        g11_val=copy(dictionary_G11[tuple_params]) # G_{1,k,q;1,k',q'}(t)

                        c00_dt=c00_val-sum_c00*dt+g00_val*dt
                        c01_dt=c01_val-sum_c01*dt+g01_val*dt
                        c10_dt=c10_val-sum_c10*dt+g10_val*dt
                        c11_dt=c11_val-sum_c11*dt+g11_val*dt
                        
                        tuple_param_dt=(t_dt,k,q,k_p,q_p) # we add the values in the next time step t+dt
                        dict_tuple_C00[tuple_param_dt]=c00_dt # C_{0,k,q;0,k',q'}(t+dt)
                        dict_tuple_C01[tuple_param_dt]=c01_dt # C_{0,k,q;1,k',q'}(t+dt)
                        dict_tuple_C10[tuple_param_dt]=c10_dt # C_{1,k,q;0,k',q'}(t+dt)
                        dict_tuple_C11[tuple_param_dt]=c11_dt # C_{1,k,q;1,k',q'}(t+dt)
                    end
                end

            end

        end

    end
    # tuples 
    open("ratio_data/C00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_tuple_C00)
    end
    open("ratio_data/C01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_tuple_C01)
    end
    open("ratio_data/C10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_tuple_C10)
    end
    open("ratio_data/C11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_tuple_C11)
    end

end

# function to compute the second derivative of β^s and β^i  over time
function SAME_second_bs_bi(type, beta,gamma,num_of_individuals, degree_val,time,dt)
    
    beta_s_array=readdlm("ratio_data/beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    beta_i_array=readdlm("ratio_data/beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]

    denominator_kq_array=readdlm("ratio_data/denominator_kq-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]
    denominator_q_array=readdlm("ratio_data/denominator_q-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt")[:,1]

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    

    #                                 t   , k , q , k', q'
    second_beta_s_tuple=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()
    second_beta_i_tuple=Dict{Tuple{Float64,Int,Int,Int,Int},Float64}()

    counter=1
    for t in 0:dt:time

        beta_s_val=beta_s_array[counter] #  β^s(t)
        beta_i_val=beta_i_array[counter] #  β^i(t)
        denominator_kq=denominator_kq_array[counter]
        denominator_q=denominator_q_array[counter]
        counter+=1

        for k in degree_array
            for q in 1:k+1
                for k_p in degree_array
                    for q_p in 1:k_p+1
                        tuple_params=(t,k,q,k_p,q_p)
                        second_beta_s_tuple[tuple_params]=second_beta_s(gamma,beta,k,q-1,k_p,q_p-1,beta_s_val,denominator_kq)
                        second_beta_i_tuple[tuple_params]=second_beta_i(gamma,beta,k,q-1,k_p,q_p-1,beta_i_val,denominator_q)
                    end
                end
            end
        end
    end
    open("ratio_data/second_beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, second_beta_s_tuple)
    end
    open("ratio_data/second_beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, second_beta_i_tuple)
    end
end

# function to compute the Hessians  over time
function SAME_PHI0_00_PHI1_00_PHI1_01_PHI1_10(type, beta,gamma,num_of_individuals, degree_val,time,dt)

    phi0_dictionary= open("ratio_data/phi0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    phi1_dictionary= open("ratio_data/phi1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    first_beta_s_dict= open("ratio_data/first_derivative_beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    first_beta_i_dict= open("ratio_data/first_derivative_beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    second_beta_s_tuple_dictionary= open("ratio_data/second_beta_s-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    second_beta_i_tuple_dictionary= open("ratio_data/second_beta_i-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    #                                t   , k , q , k', q',k'',q''
    dict_Phi0_00_tuple=Dict{Tuple{Float64,Int,Int,Int,Int,Int,Int}, Float64}()
    dict_Phi1_00_tuple=Dict{Tuple{Float64,Int,Int,Int,Int,Int,Int}, Float64}()
    dict_Phi1_01_tuple=Dict{Tuple{Float64,Int,Int,Int,Int,Int,Int}, Float64}()
    dict_Phi1_10_tuple=Dict{Tuple{Float64,Int,Int,Int,Int,Int,Int}, Float64}()


    for t in 0:dt:time
        first_beta_s_d=copy(first_beta_s_dict[t][1]) # ∂ β^s(t) for all k and q
        first_beta_i_d=copy(first_beta_i_dict[t][1]) # ∂ β^i(t) for all k and q
        for k in degree_array

            phi0_kq_array=copy(phi0_dictionary[t][1][k]) # array with ϕ_{0,k}(t) for all q
            phi1_kq_array=copy(phi1_dictionary[t][1][k]) # array with ϕ_{1,k}(t) for all q
            
            for q in 1:k+1 # k, q
                phi0=phi0_kq_array[q] # ϕ_{0,k,q}(t)
                phi1=phi1_kq_array[q] # ϕ_{1,k,q}(t)
                if q==1
                    phi0_1=0
                    phi1_1=0
                else
                    phi0_1=phi0_kq_array[q-1] # ϕ_{0,k,q-1}(t)
                    phi1_1=phi1_kq_array[q-1] # ϕ_{1,k,q-1}(t)
                end
                for k_p in degree_array # k',q'
                    for q_p in 1:k_p+1
                        first_beta_s_val_k2q2=first_beta_s_d[k_p][q_p] # ∂ β^s\∂ϕ_{0,k',q'}
                        first_beta_i_val_k2q2=first_beta_i_d[k_p][q_p] # ∂ β^i\∂ϕ_{0,k',q'}
                        for k_pp in degree_array # k'', q''
                            for q_pp in 1:k_pp+1
                                # parameters
                                tuple_params=(t,k_p,q_p,k_pp,q_pp)
                                first_beta_s_val_k3q3=first_beta_s_d[k_pp][q_pp] # ∂ β^s\∂ϕ_{0,k'',q''}
                                first_beta_i_val_k3q3=first_beta_i_d[k_pp][q_pp] # ∂ β^i\∂ϕ_{0,k'',q''}

                                second_beta_s_val=copy(second_beta_s_tuple_dictionary[tuple_params]) # ∂^2 β^s\∂ϕ_{0,k',q'}∂ϕ_{0,k'',q''}
                                second_beta_i_val=copy(second_beta_i_tuple_dictionary[tuple_params]) # ∂^2 β^i\∂ϕ_{0,k',q'}∂ϕ_{0,k'',q''}

                                # values 
                                val_Phi0_00=Phi0_00(k,k_p,k_pp,q-1,q_p-1,q_pp-1,phi0,phi0_1,first_beta_s_val_k2q2,first_beta_s_val_k3q3,second_beta_s_val)
                                val_Phi1_00=Phi1_00(k,q-1,phi1,phi1_1,second_beta_i_val)
                                val_Phi1_01=Phi1_01(k,k_pp,q-1,q_pp-1,first_beta_i_val_k2q2)
                                val_Phi1_10=Phi1_01(k,k_p,q-1,q_p-1,first_beta_i_val_k3q3)

                                # new Tuple
                                new_tuple_params=(t,k,q,k_p,q_p,k_pp,q_pp)

                                dict_Phi0_00_tuple[new_tuple_params]=val_Phi0_00
                                dict_Phi1_00_tuple[new_tuple_params]=val_Phi1_00
                                dict_Phi1_01_tuple[new_tuple_params]=val_Phi1_01
                                dict_Phi1_10_tuple[new_tuple_params]=val_Phi1_10
                            end
                        end
                    end
                end
            end
        end
    end

    open("ratio_data/Phi0_00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Phi0_00_tuple)
    end
    open("ratio_data/Phi1_00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Phi1_00_tuple)
    end
    open("ratio_data/Phi1_01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Phi1_01_tuple)
    end
    open("ratio_data/Phi1_10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Phi1_10_tuple)
    end
end

# function to compute the values of the Γ matrix over time
function SAME_GAMMA_0_GAMMA_1(type, beta,gamma,num_of_individuals, degree_val,time,dt)
    a0a0_dictionary= open("ratio_data/C00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    a0a1_dictionary= open("ratio_data/C01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    a1a0_dictionary= open("ratio_data/C10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dict_Phi0_00= open("ratio_data/Phi0_00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dict_Phi1_00= open("ratio_data/Phi1_00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dict_Phi1_01= open("ratio_data/Phi1_01-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dict_Phi1_10= open("ratio_data/Phi1_10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    
    #                    t   ,               k  ,   q=0,1,...,k
    dict_Gamma_0=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    dict_Gamma_1=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()

    for t in 0:dt:time
        dictionary_states_Gamma_0=Dict{Int64, Vector{Float64}}()
        dictionary_states_Gamma_1=Dict{Int64, Vector{Float64}}()
        
        for k in degree_array
            
            array_Gamma_0_k=[] #length: k+1
            array_Gamma_1_k=[] #length: k+1
            for q in 1:k+1
                #= sum_{k',q',k'',q''} 
                  <a_{0,k',q'}a_{0,k'',q''}> * ∂^2 Φ_{0,k,q}/(∂ϕ_{0,k',q'}∂ϕ_{0,k'',q''}) =#
                sum_Phi0=0 
                #= sum_{k',q',k'',q''} [
                  <a_{0,k',q'}a_{0,k'',q''}> * ∂^2 Φ_{1,k,q}/(∂ϕ_{0,k',q'}∂ϕ_{0,k'',q''})
                + <a_{0,k',q'}a_{1,k'',q''}> * ∂^2 Φ_{1,k,q}/(∂ϕ_{0,k',q'}∂ϕ_{1,k'',q''}) 
                + <a_{1,k',q'}a_{0,k'',q''}> * ∂^2 Φ_{1,k,q}/(∂ϕ_{1,k',q'}∂ϕ_{0,k'',q''}) ] =#
                sum_Phi1=0 
                
                for k_p in degree_array
                    
                    for q_p in 1:k_p+1
                        for k_pp in degree_array
                            for q_pp in 1:k_pp+1
                                # parameters 
                                tuple_params_aa=(t,k_p,q_p,k_pp,q_pp) # <a_{k',q'}a_{k'',q''}>
                                tuple_params_Phi=(t,k,q,k_p,q_p,k_pp,q_pp) # ∂^2Φ_{k,q}/∂ϕ_{k',q'}∂ϕ_{k'',q''}

                                a0a0=copy(a0a0_dictionary[tuple_params_aa])
                                a0a1=copy(a0a1_dictionary[tuple_params_aa])
                                a1a0=copy(a1a0_dictionary[tuple_params_aa])

                                Phi0_00_val=copy(dict_Phi0_00[tuple_params_Phi])
                                Phi1_00_val=copy(dict_Phi1_00[tuple_params_Phi])
                                Phi1_01_val=copy(dict_Phi1_01[tuple_params_Phi])
                                Phi1_10_val=copy(dict_Phi1_10[tuple_params_Phi])

                                #computations
                                sum_Phi0+=a0a0*Phi0_00_val
                                sum_Phi1+=(a0a0*Phi1_00_val+a0a1*Phi1_01_val+a1a0*Phi1_10_val)
                    
                            end
                        end
                    end
                end
                push!(array_Gamma_0_k, sum_Phi0/2)
                push!(array_Gamma_1_k, sum_Phi1/2)
            end
            dictionary_states_Gamma_0[k]=array_Gamma_0_k
            dictionary_states_Gamma_1[k]=array_Gamma_1_k
        end
        dict_Gamma_0[t]=[dictionary_states_Gamma_0]
        dict_Gamma_1[t]=[dictionary_states_Gamma_1]
    end
    open("ratio_data/Gamma_0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Gamma_0)
    end
    open("ratio_data/Gamma_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, dict_Gamma_1)
    end
end

# function to compute the values of the <b_{n,k,q}> over time
function SAME_b_equations(type, beta,gamma,num_of_individuals, degree_val,time,dt)

    dictionary_B00_tuple= open("ratio_data/B00-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_B10_tuple= open("ratio_data/B10-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_B11_tuple= open("ratio_data/B11-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_Gamma_0= open("ratio_data/Gamma_0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_Gamma_1= open("ratio_data/Gamma_1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    

    dict_b0_kq=Dict{Int64, Vector{Float64}}()
    dict_b1_kq=Dict{Int64, Vector{Float64}}()

    # we initialize a_{0,k,m}=a_{1,k,m}=0 for all k and m 

    for k in degree_array # the initial condition is all 0's
        dict_b0_kq[k]=zeros(k+1)
        dict_b1_kq[k]=zeros(k+1)
    end
    #=
    for k in degree_array # the initial condition follow a gaussian distribution
        b0_array=[]
        b1_array=[]
        for q in 1:k+1
            push!(b0_array,randn()) # mean 0 and standard deviation 1
            push!(b1_array,randn())
        end
        dict_b0_kq[k]=b0_array
        dict_b1_kq[k]=b1_array
    end
    =#

    evolution_of_b0_kq=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()
    evolution_of_b1_kq=Dict{Float64, Vector{Dict{Int64, Vector{Float64}}}}()

    evolution_of_b0_kq[0.0]=[dict_b0_kq] # the initial conditions don't follow a normal distribution
    evolution_of_b1_kq[0.0]=[dict_b1_kq]
    

    for t in 0:dt:time
        t_dt=round(t+dt, digits=number_digits)
        dictionary_states_b0=Dict{Int64, Vector{Float64}}()
        dictionary_states_b1=Dict{Int64, Vector{Float64}}()
        for k in degree_array
            b0_kq_array=copy(evolution_of_b0_kq[t][1][k]) # selects b_{0,k}(t=t_i)
            b0_updates=copy(evolution_of_b0_kq[t][1][k])

            b1_kq_array=copy(evolution_of_b1_kq[t][1][k]) 
            b1_updates=copy(evolution_of_b1_kq[t][1][k])

            Gamma_0_array=copy(dictionary_Gamma_0[t][1][k])
            Gamma_1_array=copy(dictionary_Gamma_1[t][1][k])
            for q in 1:k+1
                b0_kq=b0_kq_array[q] # b_{0,k,q}
                b1_kq=b1_kq_array[q] # b_{1,k,q}
                
                sum_b0=0
                sum_b1=0

                Gamma_0=Gamma_0_array[q] # Γ_{0,k,q} 
                Gamma_1=Gamma_1_array[q] # Γ_{1,k,q}
                
                #= loop to obtain 
                \sum_{k',q'} B_{0,k,q;0,k',q'}(t)*b_{0,k',q'}(t)
                \sum_{k',q'} B_{1,k,q;0,k',q'}(t)*b_{0,k',q'}(t)
                \sum_{k',q'} B_{1,k,q;1,k',q'}(t)*b_{1,k',q'}(t)=#
                for k_p in degree_array
                    
                    b0_kpqp_array=copy(evolution_of_b0_kq[t][1][k_p]) # selects b_{0,k'}(t=t_i)
                    b1_kpqp_array=copy(evolution_of_b1_kq[t][1][k_p]) # selects b_{1,k'}(t=t_i)
                    for q_p in 1:k_p+1
                        b0_kpqp=b0_kpqp_array[q_p] # b_{0,k',q'}
                        b1_kpqp=b1_kpqp_array[q_p] # b_{1,k',q'}
                        tuple_params=(t,k,q,k_p,q_p)
                        #                      B_{n,k,q;n',k',q'}            *b_{n',k',q'}
                        sum_b0+=dictionary_B00_tuple[tuple_params]*b0_kpqp
                        sum_b1+=dictionary_B10_tuple[tuple_params]*b0_kpqp+dictionary_B11_tuple[tuple_params]*b1_kpqp
                        #sum_B00_b0+=dictionary_B00_tuple[tuple_params]*b0_kpqp
                        #sum_B10_b0+=dictionary_B10_tuple[tuple_params]*b0_kpqp
                        #sum_B11_b1+=dictionary_B11_tuple[tuple_params]*b1_kpqp
                    end
                end

                # now we apply Euler 

                #b0_kq_dt=b0_kq-sum_B00_b0*dt+Gamma_0*dt
                #b1_kq_dt=b1_kq-sum_B10_b0*dt-sum_B11_b1*dt+Gamma_1*dt

                b0_kq_dt=b0_kq-sum_b0*dt+Gamma_0*dt
                b1_kq_dt=b1_kq-sum_b1*dt+Gamma_1*dt

                b0_updates[q]=b0_kq_dt
                b1_updates[q]=b1_kq_dt
            end

            b0_kq_array=copy(b0_updates) 
            b1_kq_array=copy(b1_updates) 

            dictionary_states_b0[k]=b0_kq_array
            dictionary_states_b1[k]=b1_kq_array

        end

        evolution_of_b0_kq[t_dt]=[dictionary_states_b0] # creates b_0 (t=t_i+dt)
        evolution_of_b1_kq[t_dt]=[dictionary_states_b1] # creates b_1 (t=t_i+dt)
    end

    open("ratio_data/b0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, evolution_of_b0_kq)
    end

    open("ratio_data/b1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "w") do io
        serialize(io, evolution_of_b1_kq)
    end
end

# function to compute the values of the <m(t)> and <ρ(t)> over time
function SAME_density_adopters_links(type, beta,gamma,num_of_individuals, degree_val,time,dt)
    dictionary_b0= open("ratio_data/b0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_b1= open("ratio_data/b1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_phi0=open("ratio_data/phi0-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end
    dictionary_phi1=open("ratio_data/phi1-$type-$num_of_individuals-$beta-$gamma-$degree_val.ser", "r") do io
        deserialize(io)
    end

    degree_array=Int.(readdlm("ratio_data/degree_array-$type-$degree_val.txt")[:,1])
    number_of_adopters_with_degree_k=readdlm("ratio_data/k_degrees-$type-$num_of_individuals-$degree_val.txt")[:,1]

    sum_degrees=0
    for k in degree_array
        position=findfirst(x->x==k, degree_array)
        sum_degrees+=k*number_of_adopters_with_degree_k[position]
    end
    
    adopters_t_array=[]
    active_links_t_array=[]
    t_array=[t for t in 0:dt:time]
    
    for t in 0:dt:time
        
        sum_phi1=0
        sum_q_phi0=0
        sum_b1=0
        sum_q_b0=0
        for k in degree_array
            
            phi0_array=copy(dictionary_phi0[t][1][k]) # ϕ_{0,k} for all q's
            phi1_array=copy(dictionary_phi1[t][1][k]) # ϕ_{1,k} for all q's
            b0_array=copy(dictionary_b0[t][1][k])     # b_{0,k} for all q's
            b1_array=copy(dictionary_b1[t][1][k])     # b_{1,k} for all q's
            for q in 1:k+1
                phi0_array[q]=(q-1)*phi0_array[q]
                b0_array[q]=(q-1)*b0_array[q]
            end

            sum_phi1+=sum(phi1_array)
            sum_q_phi0+=sum(phi0_array)
            sum_b1+=sum(b1_array)
            sum_q_b0+=sum(b0_array)
        end
        push!(adopters_t_array, sum_phi1 + sum_b1 / num_of_individuals)
        push!(active_links_t_array, num_of_individuals*(sum_q_phi0 + sum_q_b0 / num_of_individuals)*2/sum_degrees)

    end

    result_adopters=hcat(t_array, adopters_t_array)
    result_links=hcat(t_array, active_links_t_array)

    writedlm("ratio_data/SAME-adopters-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", result_adopters)
    writedlm("ratio_data/SAME-links-$type-$num_of_individuals-$beta-$gamma-$degree_val.txt", result_links)

end

function every_function_together(network,beta,gamma,Num,degree_val,time,dt)
    SAME_phi_equations(network,beta,gamma,Num,degree_val,time,dt) # check 
    SAME_first_bs_bi_and_p00_p01(network,beta,gamma,Num,degree_val,time,dt) # check
    SAME_B00_B10_B11_G00_G01_G10_G11(network,beta,gamma,Num,degree_val,time,dt)
    SAME_C00_C01_C10_C11(network,beta,gamma,Num,degree_val,time,dt)
    SAME_a_ja_k_equations(network,beta,gamma,Num,degree_val,time,dt)
    SAME_second_bs_bi(network,beta,gamma,Num,degree_val,time,dt)
    SAME_PHI0_00_PHI1_00_PHI1_01_PHI1_10(network,beta,gamma,Num,degree_val,time,dt)
    SAME_GAMMA_0_GAMMA_1(network,beta,gamma,Num,degree_val,time,dt)
    SAME_b_equations(network,beta,gamma,Num,degree_val,time,dt)
    SAME_density_adopters_links(network,beta,gamma,Num,degree_val,time,dt)
end

#= TYPES OF NETWORKS:

-square || 1
-erdos  ||0.1
-scale_free ||5
-z  ||5; 4

=#

Num=10^4

beta=1.088
gamma=0.0004
#network="square"
time=time_respect_parameters(beta, gamma, 0.999)
time_step=0.1
number_digits=count_decimals(time_step)

#every_function_together("z",beta,gamma,Num,4,time,time_step)
#every_function_together("z",beta,gamma,Num,5,time,time_step)
#every_function_together("square",beta,gamma,Num,1,time,time_step)




