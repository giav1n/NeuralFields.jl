# NeuralFields.jl

*NeuralFields.jl* is an open-source Julia based simulator of a network of interacting populations of spiking neurons.
Instead of simulating the activity of spiking neurons, we integrate for each node the membrane potential density $p(v,t)$
that follows a Fokker-Planck equation [1,2]. Inter nodes interaction is mediate via the mean and variance of the synaptic current, this allows
to implement a masivelly parallel integration.

Finite-size fluctuations are included at each node following the idesa developed in [2,3].


