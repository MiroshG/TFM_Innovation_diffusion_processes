using DifferentialEquations, Plots, LaTeXStrings, Statistics, FindPeaks1D, DelimitedFiles

# Mean Field approximation for the Bass model
function MF_approx!(du, u #a=u[1]
    , p, t)
    du[1]=(1-u[1])*(gamma+beta*u[1])
end

# Mean Field approximation for the q-Bass model
function MF_q_Bass_approx!(du, u #a=u[1]
    , p, t)
    du[1]=(1-u[1])*(gamma+beta*u[1]^q)
end

# Mean Field approximation for the NWOM-Bass model
function MF_NWOM_approx!(du, u #p=u[1] h=u[2]
    , p, t)
    du[1]=(1-u[1]-u[2])*(gamma+beta*u[1])*(1-u[2])^alpha*(1-delta)
    du[2]=(1-u[1]-u[2])*(gamma+beta*u[1])*(1-u[2])^alpha*delta
end

# Gaussian approximation for the Bass model
function gaussian_approx!(du, u #<n>=u[1], σ^2=u[2]
    , p, t)
    #du[1]=gamma*(N-u[1])+beta*(u[1]-(u[2]+u[1]^2)/N)
    #du[2]=gamma*(N-u[1]-2*u[2])+beta*(u[1]-(u[2]+u[1]^2)/N+(2-6*u[1]/N)*u[2])
    #du[2]=gamma*(N-u[1]-2*u[2])+beta*(u[1]*N+2*N*u[2]-(u[2]+u[1]^2)-4*u[1]*u[2])/N

    du[1]=gamma*(1-u[1])+beta*(u[1]-u[2]-u[1]^2)
    du[2]=gamma*(1/N-u[1]/N-2*u[2])+beta*(u[1]/N+2*u[2]-(u[2]+u[1]^2)/N-4*u[1]*u[2])
end

# function to compute the time until the density of adopters reaches a certain value a=rho for the Bass model
function time_respect_parameters(beta, gamma,rho)
    t=1/(beta+gamma)*log((1+beta*rho/gamma)/(1-rho))
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end

gamma=0.0004 # 0.3
beta=1.088 # 1.088 # 0.7 # 0.1
alpha=2 #0.01
delta=0.3
q=2
N=10^4

time=time_respect_parameters(beta, gamma, 0.999)

u0=[0.0, 0.0]
u0_MF=[0.0]
tspan_Bass=(0.0,time)
tspan_qBass=(0.0,100)
tspan_NWOMBass=(0.0,30)


function results_gaussian(gamma,beta)
    prob = ODEProblem(gaussian_approx!, u0, tspan_Bass)
    sol = solve(prob, Rosenbrock23())

    result=hcat(sol.t, sqrt.(sol[2,:]))#./N)
    writedlm("Bass_data//Gaussian_data//Gaussian_approx_std_$beta-$gamma.txt", result)

    result_2=hcat(sol.t, sol[1,:])#./N)
    writedlm("Bass_data//Gaussian_data//Gaussian_approx_mean_$beta-$gamma.txt", result_2)

    result_3=hcat(sol.t, sol[2,:])#./N^2)
    writedlm("Bass_data//Gaussian_data//Gaussian_approx_variance_$beta-$gamma.txt", result_3)

end

function results_MF(gamma,beta)
    prob = ODEProblem(MF_approx!, u0_MF, tspan_Bass)
    sol = solve(prob, Rodas5()) #Rosenbrock23())

    result=hcat(sol.t, sol[1,:])
    writedlm("Bass_data//MF_data//MF_approx_mean_$beta-$gamma.txt", result)

end

#results_gaussian(gamma,beta)
#results_MF(gamma,beta)

function results_MF_qBass(gamma,beta,q)
    prob = ODEProblem(MF_q_Bass_approx!, u0_MF, tspan_qBass)
    sol = solve(prob, Rosenbrock23()) #Rodas5()) #Rosenbrock23())

    result=hcat(sol.t, sol[1,:])
    writedlm("qBass_data//MF_qBass_approx_mean_$beta-$gamma-$q.txt", result)

    plot(sol.t, sol[1,:])
    savefig("Figures/qBass/test.png")
end

function results_MF_NWOM(gamma,beta,alpha,delta)
    prob = ODEProblem(MF_NWOM_approx!, u0, tspan_NWOMBass)
    sol = solve(prob, Rodas5()) #Rosenbrock23())

    result=hcat(sol.t, sol[1,:]+sol[2,:])
    writedlm("NWOMBass_data//MF_NWOM_adopters_$beta-$gamma-$delta-$alpha.txt", result)

    result=hcat(sol.t, sol[2,:])
    writedlm("NWOMBass_data//MF_hate_$beta-$gamma-$delta-$alpha.txt", result)

    plot(sol.t, sol[1,:]+sol[2,:])
    savefig("Figures/NWOM/test.png")

end

#results_MF_NWOM(gamma,beta,alpha,delta)
##results_MF_qBass(gamma,beta,q)
#results_Gaussian_q2_Bass(gamma,beta,q)