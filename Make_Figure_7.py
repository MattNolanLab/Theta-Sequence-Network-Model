from brian import *
from scipy import signal
from Simulate_Network_modified import Simulate_Network_modified

# Simulation Parameters

Ni = 100     # number of interneurons
Npc = 100    # number of place cells
v = 30.0 * cm / second  # running speed
f_th = 8 * Hz # theta frequency

L = 2.0 * metre

PFcent = np.zeros(Npc) 

for i in range(Npc):

	PFcent[i] = 0.5 * metre + np.random.rand(1)*1*metre # random place field mapping
	#PFcent[i] = np.random.rand(1)*5*metre 
	

# Run Network Simulation

[HPC, HI, SPC, SI, MI, MPC] = Simulate_Network_modified(Ni, Npc, v, PFcent, L, 0)


# plot results



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

subplot(211)

for i in range(Npc):

	scatter(MPC.spiketimes[i] / second, PFcent[i] * np.ones(len(MPC.spiketimes[i])))

title('Random mapping - 1 active place cell per interneuron', fontsize = 18)
ylabel('Place field center (m)', fontsize = 18)
xticks([])

axis([L / v / second / 2.0 - 0.5, L / v / second / 2.0 + 0.5, L/2.0 / metre - 0.5, L/2.0 / metre + 0.5])

# make plot with disruption

reinit(states=True)

Ni = 50

[HPC, HI, SPC, SI, MI, MPC] = Simulate_Network_modified(Ni, Npc, v, PFcent, L, 0)

subplot(212)

for i in range(Npc):

	scatter(MPC.spiketimes[i] / second, PFcent[i] * np.ones(len(MPC.spiketimes[i])))

title('Random mapping - 2 active place cells per interneuron per metre', fontsize = 18)
xlabel('Time (s)', fontsize = 18)
ylabel('Place field center (m)', fontsize = 18)

axis([L / v / second / 2.0 - 0.5, L / v / second / 2.0 + 0.5, L/2.0 / metre - 0.5, L/2.0 / metre + 0.5])

show()

