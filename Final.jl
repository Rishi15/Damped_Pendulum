using Flux, ProgressMeter
using Flux.Data: DataLoader
using Plots
using IterTools
using JLD
using Statistics
using Random
using MLDataUtils
using LinearAlgebra
using Tracker
using RDatasets
# using StatsPlots
using Flux: params
using IterTools
import Random
using ProgressBars
save_plot = false
k=300
c=5
m=1
g=9.8
batch_size = 64
include("Integrator.jl")

int=Integrator
include("Dynsys.jl")
pendulum = Dynsys.Math_pendulum(1.0, 9.81, 1.0, 300.0, 5.0, 1, 0.0)     #l, g, m, k, c, phi, phi_dot  
Random.seed!(64)
time=load("times.jld")["t"]'
time = convert(Array{Float64}, time)
Phi = load("Data.jld")["Data"]'
Phi = convert(Array{Float64}, Phi)
(Phi_train,time_train),(Phi_test,time_test)=splitobs((Phi,time); at=0.35)
    

hidden = 40
model = Chain(
    Dense(1, hidden, Flux.tanh),  
    Dense(hidden, hidden, Flux.tanh),
    Dense(hidden, hidden, Flux.tanh),
    Dense(hidden, hidden, Flux.tanh),
    Dense(hidden, 1, Flux.tanh)
)
optimizer = ADAM()

PINN=[]
dt=0.001
time_P=collect(range(0,2,step=dt))'
function res_loss(x,y)
    loss_train=Flux.mse(model(x),y)
    phi_pred = model(time_P)
    # Calculate the second derivative of phi using finite differences
    phi_dot_pred = (phi_pred[3:end] - phi_pred[1:end-2]) / (2*dt)
    phi_ddot_pred = (phi_pred[3:end] - 2*phi_pred[2:end-1] + phi_pred[1:end-2]) / (dt^2)
    r_pinn=phi_ddot_pred+5*phi_dot_pred+300*phi_pred[2:end-1]
    r_pinn=1e-4*Flux.mean(r_pinn.^2)
    return loss_train+(r_pinn)
end
trainer= DataLoader((time_train, Phi_train), shuffle=true, batchsize=40)#(t_phy,del_t)
training_loss = Float64[]
testing_loss = Float64[]
epochs = Int64[]
epochs=2000
for epochs in ProgressBar(1:epochs)
    Flux.train!(res_loss, Flux.params(model), trainer, optimizer)   # Training call
    
end
plot(time',[Phi' model(time)'],label="Training vs PINN")
if save_plot == true
    println("Saving Data Plot")
    plot(time',[Phi' model(time)'])
    savefig(case * ".png")
    println("Save Complete")
end
# anim = @animate for i ∈ 1:n
#     model(x, y, i)
# end
# gif(anim, "anim_fps15.gif", fps = 15)