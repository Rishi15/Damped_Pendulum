
using PrettyTables
using Plots, Printf
using DelimitedFiles
using Pkg
using JLD

show_video = true
println("-- pendulum euler/central diff --")
include("Dynsys.jl")
pendulum = Dynsys.Math_pendulum(1.0, 9.81, 1.0, 300.0, 5.0, 1, 0.0)    
Integ = Dynsys.Integrator(0.001, 2000)
Integ.res_phi = zeros(Integ.timesteps)
Integ.res_phi_dot = zeros(Integ.timesteps)

fig = Dynsys.create_fig(pendulum)
Dynsys.plot_state(pendulum)
display(fig)
Phi = []
times = []
for i in 1:Integ.timesteps
    fig = Dynsys.create_fig(pendulum)
    x = Dynsys.run_step(Integ, "central_diff", pendulum)
    Dynsys.plot_state(pendulum)
    display(fig)
    Integ.res_phi[i] = pendulum.phi
    Integ.res_phi_dot[i] = pendulum.phi_dot
    println(Integ.res_phi[i])
    append!(Phi, Integ.res_phi[i])            
    append!(times, i * Integ.delta_t)      
end
phaseShift = []

save("Data.jld", "Data", Phi)
save("times.jld", "t", times)

for i in 1:Integ.timesteps    
    append!(phaseShift, pi / 2)
end

plot(times,Phi)
