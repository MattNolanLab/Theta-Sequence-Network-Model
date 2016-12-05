from brian import *
from scipy import signal
from scipy import stats
from Simulate_Network import Simulate_Network
from Membrane_Theta_Freq import Membrane_Theta_Freq
from Pop_Phase_Corr import Pop_Phase_Corr
from Sequence_Corr import Sequence_Corr
from Phase_Corr import Phase_Corr


# Simulation Parameters

Npc = 1000    # number of place cells
v = 30 * cm / second  # running speed
f_th = 8 * Hz # theta frequency
Fs = 10000

L = 5.0 * metre

Ni_list = [1000, 500, 250, 200, 125, 100, 50] # different numbers of interneurons (need to be able to divide number of pyramidal cells equally by number of interneurons)

Trun = L / v

Frac_Cells = np.zeros(len(Ni_list)) # fraction of active place cells per interneuron
C = np.zeros(len(Ni_list))
IND = np.zeros(len(Ni_list))
Rho_th = np.zeros(len(Ni_list))
PF_ave = np.zeros(len(Ni_list))



# Configure place field mapping - random mapping

PFcent = np.zeros(Npc) 

for i in range(Npc):

	PFcent[i] = np.random.rand(1)*L # random place field mapping



# Simulate network and analyse phase precession while varying number of interneurons


for i in range(len(Ni_list)):

	reinit(states=True)
	
	Ni = Ni_list[i]

	# Configure place field mapping - optimal mapping
	#
	#PFcent = np.zeros(Npc) 
	#
	#for j in range(Npc/Ni):
 	#  	 for k in range(Ni):
        #
      	#		  PFcent[j + Npc/Ni*k] = j*L/(Npc/Ni) + k*L/(float(Npc/Ni))/Ni  # optimal place field mapping

	Frac_Cells[i] = float(Npc) / float(Ni)
	
	[HPC, HI, SPC, SI, MI, MPC] = Simulate_Network(Ni, Npc, v, PFcent, L)

	[C[i], IND[i]] = Pop_Phase_Corr(Npc, MPC, PFcent, v, f_th)

	Rho_th[i] = Sequence_Corr(Npc, MPC, PFcent, v, f_th, Trun, IND[i], Fs)

	PF_ave[i] = Phase_Corr(Npc, MPC, v, f_th)

	print C[i]
	print Rho_th[i]
	print PF_ave[i]



# plot results

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

plot(Frac_Cells / (L / metre), C,      color='blue',  label='Population Phase Precession',  linewidth=3.0)
plot(Frac_Cells / (L / metre), Rho_th, color='red',   label='Individual Theta Sequences',   linewidth=3.0)
plot(Frac_Cells / (L / metre), PF_ave, color='black', label='Single-Cell Phase Precession', linewidth=3.0)

legend(loc='upper right', fontsize=18)

ylabel('Correlation', fontsize=18)
xlabel('Active place cells per interneuron per meter', fontsize=18)
title('Random place field mapping')
show()





