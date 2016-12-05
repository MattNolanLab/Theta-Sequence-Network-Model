def Simulate_Network_modified( Ni, Npc, v, PFcent, L, IntNoise):

	from brian import *

    # This function uses slightly different parameters to Simulate_Network.m, and was used for the majority of simulations 

	# Note: Npc must be an integer multiple of Ni

	# Parameters of integrate and fire neuron

	tauI  = 40 * msecond      # membrane time constant, needs to be quite long to get theta spiking from dc input
	tauPC = 20 * msecond      # very variable estimates in the literature...
	Vt    = -50 * mvolt          # spike threshold
	Vr    = -70 * mvolt          # reset value (includes AHP)
	El    = -65 * mvolt          # resting potential (same as the reset)
	CmPC  = 155*pfarad          # place cell capacitance
	CmI   = 200*pfarad         # interneuron capacitance

	# Synapse time constants

	taue = 2 * msecond
	taui = 10 * msecond

	# Synapse reversal potentials

	Ee = 0 * mV
	Ei = -70 * mV

	# Other parameters

	f_th = 8*Hz   # septal theta frequency

	sigma = 40 * cm   # place field gaussian width (subthreshold field is larger than spiking field)
 
	Twidth = sigma / v  # place field width in time

	Trun = L/v # time to run along track


	# Place cell equations 

	Idc_extPC =  (1.1  + 0.005 * v / (cm / second)) * (10 ** -4)  * uamp  # input current to E cells, Gaussian amplitude 
	sigma_nE = (1.75  - 0.025* (v / (cm/second))) * mV


	eqsPC=('''
	    dVPC/dt = -(VPC - El)/tauPC - gi*(VPC - Ei)/CmPC + Idc/CmPC + sigma_nE * xi / tauPC**.5: volt
	    dgi/dt = -gi/taui        : uS
	    Idc = Idc_extPC* PF * second : amp

	    PF = exp(-(t-Tcent)**2/(2*Twidth**2)) * Hz : Hz
	    Tcent : second
	       
	    ''')


	# Interneuron equations

	Idc_MS =  0.065 * v / (cm/second) * (10 ** -6) * uamp 
	Idc0I  =  (7.95 + 0.0027 * v / (cm/second)) * (10 ** -5) * uamp 


	 

	eqsI=Equations('''
			dVI/dt = -(VI - El )/tauI - ge*(VI - Ee)/CmI + Idc/CmI + sigma_n * xi / tauI**.5: volt
	      	        dge/dt = -ge/taue        : uS
			Idc = Idc0I - Idc_MS * cos(2*pi*f_th * t)  : amp
		
			sigma_n = IntNoise * volt : volt

		       ''')
	 

	# Define cell groups

	     
	PC = NeuronGroup(N=Npc, model = eqsPC,
		  	    threshold=Vt, reset=Vr )


	for i in range(Npc):	

		PC.Tcent[i] = float(PFcent[i]) / v
	    
	   
	I = NeuronGroup(N=Ni, model = eqsI,
		    	  threshold=Vt, reset=Vr )



	# Define connections between cell groups 

	we = 0.0005*uS
	wi = 0.025*uS  # set this so as to get realistic membrane oscillations outside place field
	    
	   
	CE = Connection(PC,I,'ge')
	CI = Connection(I,PC,'gi')
	    
	# define E-I connections    
	    
	for j in range(Ni):    
	    
		CE[:,j] = 0*we/(Npc*10) # all to all component
		      
		CE[(Npc/Ni*j):(Npc/Ni*(j+1)), j] = we # phase precessing component


	# define I-E connections

	for j in range(Ni):
	    
		 CI[j, :] = 0*wi/(Ni*5) # all to all component
	
	   	 CI[j, (Npc/Ni*j):(Npc/Ni*(j+1))] = wi  # phase precessing component

	# Monitor state variables

	HPC = StateMonitor(PC, 'VPC', record=True)  
	HI = StateMonitor(I, 'VI', record=True)
	SPC = SpikeCounter(PC)
	SI = SpikeCounter(I) 
	MI = SpikeMonitor(I)
	MPC = SpikeMonitor(PC)

	   

	# Simulate
	   	

	PC.VPC = Vr # intial membrane is at resting
	I.VI = Vr   # intial membrane is at resting

	run(Trun)
	


	return (HPC, HI, SPC, SI, MI, MPC)


