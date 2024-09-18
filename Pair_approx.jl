using DifferentialEquations, Plots, LaTeXStrings, Statistics, FindPeaks1D, DelimitedFiles

function pair_approx_erdos!(du, u #σ=u[1], ρ=u[2]
    , p, t)
    mu=p[1]
    du[1]=(1-u[1])*(gamma+beta*u[2]/(2*(1-u[1])))
    du[2]=2/mu*(1-u[1])*(gamma*mu*(1-u[2]/(1-u[1])) + beta*u[2]/(1-u[1])* ( mu/2-1-(mu-1)*u[2]/(2*(1-u[1])) )  )
end

# function to compute the time until the density of adopters reaches a certain value a=rho
function time_respect_parameters(beta, gamma,rho)
    t=1/(beta+gamma)*log((1+beta*rho/gamma)/(1-rho))
    rounded_t = ceil(t)
    println(rounded_t)
    return  Int(rounded_t)
end

gamma=0.0004 # 0.0004 # 0.3
beta= 1.088 # 1-gamma # 1.088 # 0.7
N=10^4

time=time_respect_parameters(beta, gamma, 0.999)

u0=[0, 0]
tspan=(0.0,time)

function pair_approx(type, parameter)
    mu_arr=readdlm("networks_data//mean_degree_$type-$parameter-$N.txt")[1:1]
    mu=mu_arr[1]
    prob = ODEProblem(pair_approx_erdos!, u0, tspan, mu_arr)
    sol = solve(prob,Rosenbrock23())

    result=hcat(sol.t, sol[1,:])
    writedlm("ratio_data//Pair_approx-$type-$N-$beta-$gamma-$mu.txt", result)

    result_density=hcat(sol.t, sol[2,:])
    writedlm("ratio_data//Density_approx-$type-$N-$beta-$gamma-$mu.txt", result_density)
end

#=
#pair_approx("RSF", 5)
#pair_approx("erdos", 0.1)
#pair_approx("All_to_all", "-")
#pair_approx("square", 1)
#pair_approx("z", 5)
#pair_approx("z", 4)
# =#