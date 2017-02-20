from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

window_size = 1025
window_type = "hamming"
zero_padding_factor = 1
nfft = 2**(next2pow(1025) + zero_padding_factor)

parameter = Parameters(np.array([1, 3]), np.array([0.12, 0.29]), np.array([3., 0.2]))
synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, 500, 500, parameter)
synth.set_next_frame(parameter)
s = synth.get_next_frame()
w = signal.get_window(window_type, window_size)
w /= np.sum(w)
s_gt = np.zeros(window_size)
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i)*np.array(range(window_size)) +
                                             parameter.get_phase(i)*1j).real

s_gt *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = s_gt[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = s_gt[:(window_size + 1) / 2]

omega = np.arange(nfft) / float(nfft)

sfft = np.fft.fft(sw_gt, nfft)

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2,1,1)
plt.plot(omega, synth._current_spectrum._amplitude, omega, np.absolute(sfft))
plt.subplot(2,1,2)
plt.plot(omega, synth._current_spectrum._phase, omega, np.angle(sfft))
plt.figure()
plt.title("ifft of synth spectrum")
plt.plot(np.fft.ifft(synth._current_spectrum.get_complex_spectrum()))
plt.figure()
plt.title("Original and synthetized frame")
plt.plot(range(1025), s_gt, range(1025), s)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
plt.show()

plt.close()
