from brian import *
from scipy import signal
from Simulate_Network import Simulate_Network
from Membrane_Theta_Freq import Membrane_Theta_Freq

# Simulation Parameters

Ni = 1     # number of interneurons
Npc = 1    # number of place cells
v = 50 * cm / second  # running speed
f_th = 8 * Hz # theta frequency
L = 5 * metre
IntNoise = 0.15 * mV / volt 

PFcent = np.zeros(Npc) 

for i in range(Npc):

	PFcent[0] = 2.5 * metre # place field centre

# Run Network Simulation

[HPC, HI, SPC, SI, MI, MPC] = Simulate_Network(Ni, Npc, v, PFcent, L, IntNoise)


# get place field width and spike count

width = np.zeros(Npc)  # width of spiking field on this lap

for i in range(Npc):

		width[i] = float(np.amax(MPC.spiketimes[i]) - np.amin(MPC.spiketimes[i])) * v  # store average spiking field width for this session
   
print np.median(width)
print np.median(float(MPC.nspikes) / Npc)
 

# calculate spike phases and locations 

ts_PC1 = np.concatenate([MPC.spiketimes[0], MPC.spiketimes[0]])   # spike phase
ts_I   = np.concatenate([MI.spiketimes[0],   MI.spiketimes[0]])   # spike phase 

th_PC1 = np.concatenate([np.mod(MPC.spiketimes[0]*f_th*360, 360), np.mod(MPC.spiketimes[0]*f_th*360, 360) + 360])   # spike phase
th_I   = np.concatenate([np.mod(MI.spiketimes[0]*f_th*360, 360),   np.mod(MI.spiketimes[0]*f_th*360, 360) + 360])   # spike phase 
    
x_PC1 = (ts_PC1)*v  # spike location (relative to place field centre)
x_I  = ts_I*v

# Get Membrane Theta Frequencies

[F_I, Ph_I]   = Membrane_Theta_Freq(HI[0])
[F_PC, Ph_PC] = Membrane_Theta_Freq(HPC[0])

print np.amax(F_PC)
print np.amax(F_I)

# plot results



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)



subplot(423)

plot(HI.times / ms, HI[0] / mV, color = 'red', linewidth=2.0)
#ylabel('Membrane Potential (mV)', fontsize = 18)
#title('Interneuron Membrane Trace', fontsize = 22)
axis([(PFcent / v)/ ms - 2 * 1000, (PFcent / v) / ms + 2 * 1000, -85, 0])
xticks([])
yticks([])

T = MI.spiketimes[0] / ms

for i in range(len(T)):

	plot([T[i], T[i]], [-45, - 10], color = 'red', linewidth=2.0)

subplot(422)

plot(HI.times / ms, F_I, color = 'red', linewidth=2.0)

ylabel('Freq (Hz)', fontsize = 18)
axis([(PFcent / v)/ ms - 2 * 1000, (PFcent / v) / ms + 2 * 1000, 7.5, 10])
xticks([])

#plot(HPC.times / ms, - 40 - 2 * cos(2*pi*f_th * HPC.times), color='black') # theta pacemaker

subplot(427)

T = MPC.spiketimes[0] / ms

for i in range(len(T)):

	plot([T[i], T[i]], [-45, - 10], color = 'blue', linewidth=2.0)


plot(HPC.times / ms, HPC[0] / mV, color = 'blue', linewidth=2.0)

xlabel('Time (ms)', fontsize = 18)
#ylabel('Membrane Potential (mV)', fontsize = 18)
#title('Place Cell Membrane Trace',fontsize = 22)
axis([(PFcent / v)/ ms - 2 * 1000, (PFcent / v) / ms + 2 * 1000, -85, 0])
yticks([])

subplot(426)

plot(HPC.times / ms, F_PC, color = 'blue', linewidth=2.0)

ylabel('Freq (Hz)', fontsize = 18)
axis([(PFcent / v)/ ms - 2 * 1000, (PFcent / v) / ms + 2 * 1000, 7.5, 10])
xticks([])

subplot(424)

scatter(x_I/cm, th_I - 360, color = 'red')
ylabel('Spike Phase', fontsize = 18)
#title('Interneuron Spike Phase', fontsize = 22)
axis([150, 350, -400, 400])
yticks([-360, -180, 0, 180, 360])
xticks([])

subplot(428)

scatter(x_PC1/cm, th_PC1 - 360, color = 'blue')
xlabel('Location (cm)', fontsize = 18)
ylabel('Spike Phase', fontsize = 18)
#title('Place Cell Spike Phase', fontsize = 22)
axis([150, 350, -400, 400])
yticks([-360, -180, 0, 180, 360])


subplot(421)

plot(HPC.times / ms, cos(2*pi*f_th * HPC.times), color='black', linewidth=2.0) # theta pacemaker
axis([(PFcent / v)/ ms - 2 * 1000, (PFcent / v) / ms + 2 * 1000, -5, 5])
xticks([])
yticks([])


show()

