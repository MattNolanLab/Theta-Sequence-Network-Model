def Membrane_Theta_Freq(Trace):

	from scipy import signal
	from brian import *

	"calculate instantaneous phase and frequency of a time series"
	
	Fs = 10000 # sampling rate (Hz)
	b, a = signal.butter(2, [0.00125, 0.002], btype='band', analog=False)
	Nsmooth = 5000 # smoothing kernel width in pixels

	S = signal.filtfilt(b, a, Trace)
		
	H = signal.hilbert(S)  # hilbert transform of membrane potential
	Phase0 = np.angle(H)

	Phase_smooth = np.convolve(np.unwrap(Phase0), np.ones(Nsmooth)/Nsmooth, 'same')
	freq_smooth = np.gradient(Phase_smooth, 1) * np.divide(Fs, 2*np.pi)

	return (freq_smooth, Phase_smooth)

