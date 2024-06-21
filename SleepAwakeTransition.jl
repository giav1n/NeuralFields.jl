using Plots,ProgressMeter,DelimitedFiles
include("NeuralFields.jl")

## Sleep-Awake transition


NV=200 # Number of grid ponints 
α=0.5  # Vmin=-α*θ
dt=0.001 #[s]
cm,S,μx,σ2x,g,Npop,KJ,KJ2=NeuralFields.InitializeFPfromPerseusNet("D3.6/m_cortsurf.ini","D3.6/c_cortsurf.ini",NV,dt,α);

##
Life=25 #[s]
steps=Int(round(Life/dt))
ν,νt,c=zeros(Npop,steps),zeros(Npop),zeros(Npop)
ϵt=0.8
T0=20

function ϵ(t)
    if t<T0
        return (1-ϵt)+(ϵt)*(T0-t)/T0
    else
        return 1- ϵt
    end
end    


@showprogress for n=2:steps
    gt=g.*ϵ(n*dt)
    # Compute the recurrent contribution
    μ= KJ*νt + μx -gt.*c  # Mean of synaptic currents
    σ2=KJ2*νt +σ2x       # Variance of synaptic currents
    # c is the adaptation variable dc/dt =-c(t)/τC +ν(t)

    Threads.@threads for i=1:Npop
        NeuralFields.IntegrateFP!(cm[i],μ[i], σ2[i],S[i])
        ν[i,n],νt[i],c[i]=S[i][cm[i].NV+1:cm[i].NV+3]
   end
   

end


##

p1=heatmap(log.(ν[1:3:Npop,1:200:end] .+1),clims=(log(1),log(101)),xticks=[],ylabel=" # id of neuron",colorbar=false)
t=range(0,Life,length=steps)
p2=plot(t[1:200:end],g[1].*ϵ.(t[1:200:end]),color=:black,label=false,ylim=(0,11),yticks=[0,5,10])
plot(p1,p2, layout=grid(2,1, heights=(6/8,2/8)))
savefig("AwakeTransit.svg")