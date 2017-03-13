from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt
from utils import *

fs = 44100
window_length = 23e-3

analysis_hop, synthesis_hop = 505, 505

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 4
nfft = 2**(next2pow(window_size) + zero_padding_factor)

parameter = NonStationaryParameters(amplitudes=np.array([1.]), frequencies=np.array([0.1]), phases=np.array([3.]),
                                    acrs=np.array([0.]), fcrs=np.array([8000.]))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

n = np.arange(int(round(window_size + analysis_hop)))
t = np.arange(int(round(window_size + analysis_hop)))/float(fs)

s_gt = np.zeros(int(round(window_size + analysis_hop)))
for i in xrange(parameter.get_number_sinuses()):
    s_gt += (parameter.get_amplitude(i) + parameter.get_acrs(i)*t)*np.exp(1j*(2*np.pi*parameter.get_frequency(i) *
                                                                          n + parameter.get_phase(i) +
                                                                          parameter.get_fcrs(i) * t**2 / 2.)).real

frame = np.zeros(window_size)
frame[:window_size] = s_gt[:window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

sfft = np.fft.fft(sw_gt, nfft)
phases = np.array([np.angle(sfft[int(round(parameter.get_frequency(k) * nfft))]) -
                   (2*np.pi*parameter.get_frequency(k)*analysis_hop) +
                   (parameter.get_fcrs(k) * analysis_hop**2 / 2.)
                   for k in xrange(parameter.get_number_sinuses())])
parameter._phases = phases
synth = NonStationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter,
                                 True, (-1e-5, 1e-5), (-1e5, 1e5), 20, 20)
synth.set_next_frame(parameter)
s = synth.get_next_frame()

# Comparison
omega = np.arange(nfft) / float(nfft)

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2, 1, 1)
plt.xlabel('Normalized frequencies')
plt.plot(omega, 10*np.log10(synth.current_spectrum.amplitude**2), omega, 10*np.log10(np.absolute(sfft)**2), "r--")
plt.legend(['Synthesised spectrum', 'Actual spectrum'])
plt.title("Synthetized and real spectrums")
plt.subplot(2, 1, 2)
plt.plot(omega, np.mod(synth.current_spectrum.phase, 2*np.pi) - np.pi, omega, np.angle(sfft))
plt.figure()
plt.title(["Synthetised and original frame"])
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

parameter._phases += np.array([(2 * np.pi * parameter.get_frequency(k) * analysis_hop) +
                                   (parameter.get_fcrs(k) * analysis_hop ** 2 / 2.)
                                   for k in xrange(parameter.get_number_sinuses())])
parameter._frequencies += np.array([(analysis_hop * a) / (2 * np.pi * fs ** 2) for a in parameter._fcrs])
parameter._amplitudes += np.array([(analysis_hop * a) / float(fs) for a in parameter._acrs])
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
