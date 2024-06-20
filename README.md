# NeuralFields.jl
## The general idea
*NeuralFields.jl* is an open-source Julia based simulator of a network of interacting populations of spiking neurons.
Instead of simulating the activity of spiking neurons, we integrate for each node the membrane potential density $p(v,t)$
that follows a Fokker-Planck equation [1,2]. Inter nodes interaction is mediate via the mean and variance of the synaptic current, this allows
to implement a masivelly parallel integration.
Finite-size fluctuations are included at each node following the idesa developed in [2,3].

For the networks parameter and synaptic connections we use the initialization file Modules.ini and Connectivity.ini of the PERSEO spiking network simulator [4]

If you used this tool in your work consider to cite us in [3]

## References
[1] Brunel, Nicolas, and Vincent Hakim. "Fast global oscillations in networks of integrate-and-fire neurons with low firing rates." Neural computation 11.7 (1999): 1621-1671.
[2] Mattia, Maurizio, and Paolo Del Giudice. "Population dynamics of interacting spiking neurons." Physical Review E 66.5 (2002): 051917.
[3] Vinci, Gianni V., Roberto Benzi, and Maurizio Mattia. "Self-consistent stochastic dynamics for finite-size networks of spiking neurons." Physical Review Letters 130.9 (2023): 097402.
[4] Mattia, Maurizio, and Paolo Del Giudice. "Efficient event-driven simulation of large networks of spiking neurons and dynamical synapses." Neural computation 12.10 (2000): 2305-2329.

## How to use
All the required functions are contained in the module NeuralField.jl. An exemple of a large neural field with step by stepe guidance can be found in Example.jl 



