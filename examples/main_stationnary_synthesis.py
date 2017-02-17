from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

parameter = Parameters(np.array([1]), np.array([0.12]), np.array([1]))
synth = StationarySynthesizer(1015, "hamming", 3, 500, 500, parameter)
synth.set_next_frame(parameter)
s = synth.get_next_frame()  # todo : Issues withe the phase vocoder
w = signal.get_window("hamming", 1015)
s_gt = np.array([1])*np.exp(2*np.pi*1j*0.12*np.array(range(1015)) + 1j)
s_gt *= w
plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2,1,1)
plt.plot(synth._current_spectrum._amplitude)
plt.subplot(2,1,2)
plt.plot(synth._current_spectrum._phase)
plt.figure()
plt.title("ifft of synth spectrum")
plt.plot(np.fft.ifft(synth._current_spectrum.get_complex_spectrum()))
plt.figure()
plt.title("Original and synthetized frame")
plt.plot(range(1015), s_gt, range(1015), s)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
plt.show()

plt.close()
