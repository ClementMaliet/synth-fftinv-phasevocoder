from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
window_length = 23e-3

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 0
nfft = 2**(next2pow(window_size) + zero_padding_factor)

phases = []
phases_gt = []
frequencies = np.linspace(0, 0.5, 1000)

parameter = Parameters(np.array([1]), np.array([0.16]), np.array([1.2]))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, 100, 100, parameter)

for freq in frequencies:
    parameter._frequencies[0] = freq
    s_gt = np.zeros(window_size)
    for i in xrange(parameter.get_number_sinuses()):
        s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i)*np.array(range(window_size)) +
                                                  parameter.get_phase(i)*1j).real

    s_gt *= w

    sw_gt = np.zeros(nfft)
    sw_gt[:(window_size - 1) / 2] = s_gt[((window_size + 1) / 2):]
    sw_gt[nfft - (window_size - 1) / 2 - 1:] = s_gt[:(window_size + 1) / 2]

    sfft = np.fft.fft(sw_gt, nfft)
    phase = np.unwrap(np.angle(sfft))
    phases_gt.append(np.array([phase[int(round(a * nfft))] for a in parameter._frequencies]))

    # Synthesis
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()
    phases.append(np.array([synth._current_spectrum._phase[int(round(a * nfft))] for a in parameter._frequencies]))

# Plots
plt.figure()
plt.plot(frequencies, phases_gt, frequencies, phases)
plt.xlabel("Normalised frequency")
plt.ylabel("Sinus phase")


plt.figure()
plt.plot(frequencies, phases_gt)
plt.xlabel("Normalised frequency")
plt.ylabel("Sinus phase")
plt.show()
