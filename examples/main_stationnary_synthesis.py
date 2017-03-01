from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
window_length = 23e-3

analysis_hop, synthesis_hop = 505, 505

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 2
nfft = 2**(next2pow(window_size) + zero_padding_factor)

parameter = Parameters(np.array([500.]), np.array([0.05]), np.array([3.]))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

s_gt = np.zeros(window_size + analysis_hop)
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i) *
                                              np.array(range(window_size+analysis_hop)) +
                                              parameter.get_phase(i)*1j).real

frame = np.zeros(window_size)
frame[:window_size] = s_gt[:window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

sfft = np.fft.fft(sw_gt, nfft)
phases = np.array([np.angle(sfft[int(round(a * nfft))]) - 2*np.pi*a*analysis_hop for a in parameter._frequencies])


# Synthesis
parameter._phases = phases
synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)
synth.set_next_frame(parameter)
s = synth.get_next_frame()

# Comparison
omega = np.arange(nfft) / float(nfft)

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2, 1, 1)
plt.plot(omega, 10*np.log10(synth.current_spectrum.amplitude**2), omega, 10*np.log10(np.absolute(sfft)**2))
plt.subplot(2, 1, 2)
plt.plot(omega, np.mod(synth.current_spectrum.phase, 2*np.pi) - np.pi, omega, np.angle(sfft))
plt.figure()
plt.title("Synthetised and original frame")
plt.plot(range(window_size), s, range(window_size), frame)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
print "RMSE = " + str((1./window_size)*np.sqrt(np.sum((s-frame)**2)))

# Second frame
# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)


frame[:window_size] = s_gt[analysis_hop:analysis_hop + window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

sfft = np.fft.fft(sw_gt, nfft)

phases = np.array([synth.current_spectrum.phase[int(round(a * nfft))] for a in parameter._frequencies])

# Synthesis
parameter._phases = phases
synth.set_next_frame(parameter)
s = synth.get_next_frame()

# Comparison

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2, 1, 1)
plt.plot(omega, 10*np.log10(synth.current_spectrum.amplitude**2), omega, 10*np.log10(np.absolute(sfft)**2))
plt.subplot(2, 1, 2)
plt.plot(omega, np.mod(synth.current_spectrum.phase, 2*np.pi) - np.pi, omega, np.angle(sfft))
plt.figure()
plt.title("Synthetised and original frame")
plt.plot(range(window_size), s, range(window_size), frame)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
print "RMSE = " + str((1./window_size)*np.sqrt(np.sum((s-frame)**2)))

plt.show()
