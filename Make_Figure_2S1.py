
from brian import *

# This script simulates the reduced oscillator model, with a slowly varying square-pulse, upward ramp, or downward ramp input drive. 
# Varying current amplitudes are simulated and the resulting phase precession curves are plotted.

dt = 0.1 * ms # timestep

T = 4.0 * second # simulation time

Nstep = np.int(T / dt) # number of time steps

t = np.linspace(0, T / second, Nstep)

v = 40 * cm / second # running speed

sigma = 10 * cm # with of Gaussian drive
Tbox = 0.5 * second # width of square-pulse drive

omega_pace = 2 * pi * 8 * Hz # pacemaker frequency
Asynch = 2 * pi * Hz # synchronisation factor

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

# set colourmap 

parameters = np.linspace(0,25,26)

norm = matplotlib.colors.Normalize(
    vmin=np.min(parameters),
    vmax=np.max(parameters))

c_m = matplotlib.cm.spring

s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])


# numerically simulate differential equation

for j in range(25):
        
## Square pulse drive
        
        domega = 2 * pi * (j+1) * 0.075 * Hz # max detuning (i.e., frequency difference between pacemaker and interneuron at centre of input field). For a linear f-I curve, this is proportional to current input.   
        omega_free = omega_pace + domega * (abs(t-T/2) < Tbox) # free running oscillator frequency
        phi = np.zeros(Nstep) # oscillator phase
        phi[0] = 0

        for i in range(Nstep-1):
    
            phi[i+1] = phi[i] + (omega_pace - omega_free[i] * Hz) * dt  - Asynch * np.sin(phi[i]) * dt
        
        subplot(231)
        plot(t / second, omega_free, linewidth=2.0, color=s_m.to_rgba(j))
        subplot(234)
        plot(t / second, phi * 180 / pi, linewidth=2.0, color=s_m.to_rgba(j))
       
## ramp up drive

        domega = 2 * pi * (j+1) * 0.075 * Hz        
        omega_free = omega_pace + domega * (abs(t-T/2) < Tbox) * (t - (T/2 - Tbox)) / Tbox # free running oscillator frequency
        phi = np.zeros(Nstep) # oscillator phase
        phi[0] = 0

        for i in range(Nstep-1):
    
            phi[i+1] = phi[i] + (omega_pace - omega_free[i] * Hz) * dt  - Asynch * np.sin(phi[i]) * dt
        
        subplot(232)
        plot(t / second, omega_free, linewidth=2.0, color=s_m.to_rgba(j))
        subplot(235)
        plot(t / second, phi * 180 / pi, linewidth=2.0, color=s_m.to_rgba(j))


## ramp down drive

        domega = 2 * pi * (j+1) * 0.075 * Hz        
        omega_free = omega_pace + domega * (abs(t-T/2) < Tbox) * ((T/2 - Tbox) - t + 2*Tbox) / Tbox # free running oscillator frequency
        phi = np.zeros(Nstep) # oscillator phase
        phi[0] = 0

        for i in range(Nstep-1):
    
            phi[i+1] = phi[i] + (omega_pace - omega_free[i] * Hz) * dt  - Asynch * np.sin(phi[i]) * dt
        
        subplot(233)
        plot(t / second, omega_free, linewidth=2.0, color=s_m.to_rgba(j))
        subplot(236)
        plot(t / second, phi * 180 / pi, linewidth=2.0, color=s_m.to_rgba(j))

        
subplot(231)
xticks([])
yticks([])
ylabel('Input Current')     

axis([0, T / second, 2*pi*8 - pi, 2*pi*8 + 8*pi])   


subplot(232)
xticks([])
yticks([])
ylabel('Input Current')     

axis([0, T / second, 2*pi*8 - pi, 2*pi*8 + 8*pi])

subplot(233)
xticks([])
yticks([])
ylabel('Input Current')     

axis([0, T / second, 2*pi*8 - pi, 2*pi*8 + 8*pi])

subplot(234)        
ylabel('Phase (degrees)', fontsize = 18)
xlabel('Time (seconds)', fontsize = 18)
axis([0, T / second, -550, 100])
yticks([-360, -180, 0])

subplot(235)        
ylabel('Phase (degrees)', fontsize = 18)
xlabel('Time (seconds)', fontsize = 18)
axis([0, T / second, -550, 100])
yticks([-360, -180, 0])

subplot(236)        
ylabel('Phase (degrees)', fontsize = 18)
xlabel('Time (seconds)', fontsize = 18)
axis([0, T / second, -550, 100])
yticks([-360, -180, 0])
