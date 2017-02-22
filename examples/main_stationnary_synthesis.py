from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
window_length = 23e-3

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hamming"
zero_padding_factor = 3
nfft = 2**(next2pow(window_size) + zero_padding_factor)

parameter = Parameters(np.array([1]), np.array([0.05]), np.array([3.]))
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

omega = np.arange(nfft) / float(nfft) - 0.5

sfft = np.fft.fftshift(np.fft.fft(sw_gt, nfft))

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2,1,1)
plt.plot(omega, 10*np.log10(np.fft.fftshift(synth._current_spectrum._amplitude)**2), omega, 10*np.log10(np.absolute(sfft)**2))
plt.subplot(2,1,2)
plt.plot(omega, np.fft.fftshift(synth._current_spectrum._phase), omega, 2*np.angle(sfft))
plt.figure()
plt.title("Synthetised and original frame")
plt.plot(range(window_size), s, range(window_size), s_gt, )
plt.xlabel("Echantillons")
plt.ylabel("Signal")
plt.show()

plt.close()
