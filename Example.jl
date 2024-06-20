using Plots,ProgressMeter
include("NeuralFields.jl")
#Load parameters from PERSEO ini files
NV=200 # Number of grid ponints 
α=0.5  # Vmin=-α*θ
dt=0.001 #[s]
cm,S,μx,σ2x,g,Npop,KJ,KJ2=NeuralFields.InitializeFPfromPerseusNet("m_cortsurf.ini","c_cortsurf.ini",NV,dt,α);

##
Life=10 #[s]
steps=Int(round(Life/dt))
ν,νt,c=zeros(Npop,steps),zeros(Npop),zeros(Npop)

@showprogress for n=2:steps
    # Compute the recurrent contribution
    μ= KJ*νt + μx -g.*c  # Mean of synaptic currents
    σ2=KJ2*νt +σ2x       # Variance of synaptic currents
    # c is the adaptation variable dc/dt =-c(t)/τC +ν(t)

    Threads.@threads for i=1:Npop
        NeuralFields.IntegrateFP!(cm[i],μ[i], σ2[i],S[i])
        ν[i,n],νt[i],c[i]=S[i][cm[i].NV+1:cm[i].NV+3]
   end
   

end

heatmap(log.(ν[1:3:Npop,1:200:end] .+1),clims=(log(1),log(101)),xticks=[],xlabel="time [s]",ylabel=" # id of neuron")


