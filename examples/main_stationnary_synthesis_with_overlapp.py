#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
sine_duration = 100e-3  # s
window_length = 23e-3  # s

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 0
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 305
synthesis_hop = 305

# parameter = Parameters(np.array([1, 0.2, 2, 1.5]), np.array([0.1, 0.3, 0.23, 0.11]), np.array([3., -2., 5, -0.3]))
parameter = Parameters(np.array([1.]), np.array([0.05]), np.array([3.]))

# Find the max number of slices that can be obtained
number_frames = int(np.floor((sine_duration*fs-window_size)/analysis_hop))

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)

s_gt = np.zeros(int(round(sine_duration*fs)))
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i) *
                                              np.arange(int(round(sine_duration*fs))) +
                                              parameter.get_phase(i)*1j).real

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

    phases = np.array([synth.current_spectrum.phase[int(round(a * nfft))] for a in parameter._frequencies])
    parameter._phases = phases  # We feed the phase advance back to the synthesizer


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
                                                            s_gt[:min(len(s_gt), len(vector_time))])**2)))

plt.show()



'''
##################################
IN MAIN :

%% Finalize

% Overlap add in a vector
outputTimeStretched = fusionFrames(outputy,hopOut);

% Resample with linearinterpolation
outputTime = interp1((0:(length(outputTimeStretched)-1)),outputTimeStretched,(0:alpha:(length(outputTimeStretched)-1)),'linear');

% Return the result
outputVector = outputTime;

return

#################################
IN FUSION FRAMES FUNCTION :
framesMatrix  :  Matrix made of all the frames

function vectorTime = fusionFrames(framesMatrix, synthesis_hop)

sizeMatrix = size(framesMatrix);

% Get the number of frames
numberFrames = sizeMatrix(1);

% Get the size of each frame
sizeFrames = sizeMatrix(2);

% Define an empty vector to receive result
vectorTime = zeros(numberFrames*synthesis_hop-synthesis_hop+sizeFrames,1);

timeIndex = 1;

% Loop for each frame and overlap-add
for index=1:numberFrames

    vectorTime(timeIndex:timeIndex+sizeFrames-1) = vectorTime(timeIndex:timeIndex+sizeFrames-1) + framesMatrix(index,:)';

    timeIndex = timeIndex + synthesis_hop;

end

'''