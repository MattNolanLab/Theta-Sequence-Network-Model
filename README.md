# Theta-Sequence-Network-Model
Code for the 2016 eLife paper

Simulates the basic circuit model to produce Figure 1, the reduced oscillator model to produce Figure 2 and Figure 2-figure Supplement 1, and the full network model to produce Figure 7.

File descriptions:

Simulate_Network - this function simulates the basic circuit shown in Figure 1. It is called from Make_Figure_1

Simulate_Network_modified - this is a slightly modified version of Simulate_network, with some parameter changes.

Membrane_Theta_Freq - this function estimates the instantaneous theta frequency of an input signal. It is called from Make_Figure_1

Make_Figure_1 - simulates and plots the results of Figure 1 from the paper (not that noise varies on each simulation)

Make_Figure_2 - simulates and plots the results of Figure 2

Make_Figure_2S1 - simulates and plots the results of Figure 2-figure supplement 1

Make_Figure_7 - simulates part A from Figure 7

Make_Figure_7_sweep - sweeps through networks with optimal or random maps of different place field densities, and plots the theta sequence and phase precession statistics under each map. This needs to be run multiple times and averaged to obtain a reliable average as shown in the paper.
