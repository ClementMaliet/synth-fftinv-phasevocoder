#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt
from utils import gen_triangle_parameters

fs = 44100
sine_duration = 96e-3  # s
window_length = 23e-3  # s

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 0
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 205
synthesis_hop = 205

parameter = NonStationaryParameters(np.array([1.]), np.array([0.05]), np.array([3.]), np.array([-0.5]), np.array([1000]))

# Find the max number of slices that can be obtained
number_frames = int(np.floor((sine_duration*fs-window_size)/analysis_hop))

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)

n = np.arange(int(round(sine_duration*fs)))
t = np.arange(int(round(sine_duration*fs)))/float(fs)

s_gt = np.zeros(int(round(sine_duration*fs)))
for i in xrange(parameter.get_number_sinuses()):
    s_gt += (parameter.get_amplitude(i) + parameter.get_acrs(i)*t)*np.exp(1j*(2*np.pi*parameter.get_frequency(i) *
                                                                          n + parameter.get_phase(i) +
                                                                          parameter.get_fcrs(i) * t**2 / 2)).real

# We compute the true initial phase to take into account the windowing effect on phase
frame = np.zeros(window_size)
frame[:window_size] = s_gt[:window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

sfft = np.fft.fft(sw_gt, nfft)
phases = np.array([np.angle(sfft[int(round(a * nfft))]) - 2*np.pi*a*analysis_hop for a in parameter._frequencies])
parameter._phases = phases

for i in xrange(number_frames):
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()
    vector_frames[i, :] = s
    parameter._phases += np.array([2 * np.pi * analysis_hop *
                                   (parameter.get_frequency(k) + (analysis_hop * parameter.get_fcrs(k)) /
                                    (4 * np.pi * fs**2)) - (2*np.pi*parameter.get_frequency(k)*analysis_hop)
                                   for k in xrange(parameter.get_number_sinuses())])
    parameter._frequencies += np.array([(analysis_hop * a) / (4 * np.pi * fs ** 2) for a in parameter._fcrs])
    parameter._amplitudes += np.array([(analysis_hop * a) / float(fs) for a in parameter._acrs])
    # phases = np.array([synth.current_spectrum.phase[int(round(a * nfft))] for a in parameter._frequencies])
    # parameter._phases = phases  # We feed the phase advance back to the synthesizer


# Define an empty vector to receive result
vector_time = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)
s_gt_env = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)

# Loop for each frame and overlap-add
time_index = 0
for i in xrange(number_frames):
    vector_time[time_index:time_index + window_size] += vector_frames[i, :]*np.sum(w)
    s_gt_env[time_index:time_index + window_size] += w*np.sum(w)
    time_index += synthesis_hop

# Comparison
plt.title("Synthetised and original frame")
plt.plot(range(len(vector_time)), vector_time, label="Synthesized")
plt.plot(range(min(len(s_gt), len(vector_time))), s_gt_env[:min(len(s_gt), len(vector_time))]*s_gt[:min(len(s_gt), len(vector_time))], '--', label="Original")
plt.xlabel("Echantillons")
plt.ylabel("Signal")
plt.legend()
print "RMSE = " + str((1./len(vector_time))*np.sqrt(np.sum((vector_time[:min(len(s_gt), len(vector_time))] -
                                                            s_gt_env*s_gt[:min(len(s_gt), len(vector_time))])**2)))

plt.show()
