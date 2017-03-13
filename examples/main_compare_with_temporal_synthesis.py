#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt
from utils import *
import time

fs = 44100
sine_duration = 1  # s
window_length = 23e-3  # s

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 0
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 505
synthesis_hop = 505

amp = 3.
freq = 0.05
phi = 3.
n_harmonic = 10

parameter = gen_triangle_parameters(amp, freq, phi, n_harmonic)

# Find the max number of slices that can be obtained
number_frames = int(np.floor((sine_duration*fs-window_size)/analysis_hop))

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)

s_gt = np.zeros(window_size)
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i) *
                                              np.arange(window_size) +
                                              parameter.get_phase(i)*1j).real

# We compute the true initial phase to take into account the windowing effect on phase
frame = np.zeros(window_size)
frame[:window_size] = s_gt[:window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

sfft = np.fft.fft(sw_gt, nfft)


# FTT-1 Synthesis
phases = np.array([np.angle(sfft[int(round(a * nfft))]) - 2*np.pi*a*analysis_hop for a in parameter._frequencies])
parameter._phases = phases

print "FFT-1 synthesis"

t0 = time.clock()
for i in xrange(number_frames):
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()
    vector_frames[i, :] = s
    parameter._phases += np.array([2 * np.pi * analysis_hop * a for a in parameter._frequencies])
    # We feed the phase advance back to the synthesizer


# Define an empty vector to receive result
vector_time = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)

# Loop for each frame and overlap-add
time_index = 0
for i in xrange(number_frames):
    vector_time[time_index:time_index + window_size] += vector_frames[i, :]*np.sum(w)
    time_index += synthesis_hop
t0 = time.clock() - t0


# Time-domain synthesis

parameter = gen_triangle_parameters(amp, freq, phi, n_harmonic)
vector_frames_t = np.zeros((number_frames, window_size))

print "Time domain synthesis"

t1 = time.clock()
for i in xrange(number_frames):
    s = np.zeros(window_size)
    for j in xrange(parameter.get_number_sinuses()):
        s += parameter.get_amplitude(j) * np.exp(2 * np.pi * 1j * parameter.get_frequency(j) *
                                                 np.arange(window_size) +
                                                 parameter.get_phase(j) * 1j).real
    vector_frames_t[i, :] = s*w
    parameter._phases += np.array([2*np.pi*analysis_hop*a for a in parameter._frequencies])
    # We feed the phase advance back to the synthesizer


# Define an empty vector to receive result
vector_time_t = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)

# Loop for each frame and overlap-add
time_index = 0
for i in xrange(number_frames):
    vector_time_t[time_index:time_index + window_size] += vector_frames_t[i, :]*np.sum(w)
    time_index += synthesis_hop
t1 = time.clock() - t1

# Additive synthesis
parameter = gen_triangle_parameters(amp, freq, phi, n_harmonic)

print "Reference synthesis"

t2 = time.clock()
s_gt = np.zeros(int(round(sine_duration*fs)))
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i) *
                                              np.arange(int(round(sine_duration*fs))) +
                                              parameter.get_phase(i)*1j).real
t2 = time.clock() - t2

s_gt_env = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)
# Loop for each frame and overlap-add
time_index = 0
for i in xrange(number_frames):
    s_gt_env[time_index:time_index + window_size] += w*np.sum(w)

print "----------------------|| Benchmark ||------------------------"
print "We synthesised a signal " + str(sine_duration) + "s long, comprising in " + str(number_frames) + " frames."
print "The signal comprised in " + str(n_harmonic) + " harmonics."
print "FFT-1 synthesis : " + str(t0) + "s"
print "Time domain synthesis (frames) : " + str(t1) + "s"
print "Time domain synthesis (additive synthesis) : " + str(t2) + "s"


print "RMSEfft-1 = " + str((1./len(vector_time))*np.sqrt(np.sum((vector_time[:min(len(s_gt), len(vector_time))] -
                                                            s_gt_env * s_gt[:min(len(s_gt), len(vector_time))])**2)))

print "RMSEtemporal = " + str((1./len(vector_time_t))*np.sqrt(np.sum((vector_time_t[:min(len(s_gt), len(vector_time_t))] -
                                                            s_gt_env * s_gt[:min(len(s_gt), len(vector_time_t))])**2)))