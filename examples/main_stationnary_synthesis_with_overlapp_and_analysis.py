from sources.Synthesizer import *
import matplotlib.pyplot as plt
from estimation_func import *
from numpy import *
from matplotlib.pyplot import *
from scipy import signal
from utils import *
from scipy.io import wavfile

filename = "Mix.wav"
window_length = 23e-3  # s

[fs, s_gt] = wavfile.read(filename)
sine_duration = 46e-3#len(s_gt)/fs  # s

parameter = Parameters(np.array([500.]), np.array([0.05]), np.array([3.]))
s_gt = np.zeros(int(round(sine_duration*fs)))
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i) *
                                              np.arange(int(round(sine_duration*fs))) +
                                              parameter.get_phase(i)*1j).real

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 2
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 505
synthesis_hop = 505

w = signal.get_window(window_type, window_size)
w /= np.sum(w)

params = ParametersBuilder(s_gt, fs, analysis_hop, window_length, window_type, zero_padding_factor, 40)

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                              Parameters.void_parameters(1))
# Find the max number of slices that can be obtained
number_frames = params.number_frames

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# # Analysis
# w = signal.get_window(window_type, window_size)
# w /= np.sum(w)
#
# # We compute the true initial phase to take into account the windowing effect on phase
# frame = np.zeros(window_size)
# frame[:window_size] = s_gt[:window_size]
# frame *= w
#
# sw_gt = np.zeros(nfft)
# sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
# sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]
#
# [a1, omega1, phi1, acr, fcr] = estimator1(frame, fs, w, nfft, th)
# freq1 = omega1 /(2*np.pi*fs)
# parameter = Parameters(a1, freq1, phi1)

for i in xrange(number_frames):
    # frame[:window_size] = s_gt[i*analysis_hop:i*analysis_hop + window_size]
    #
    # [a1, omega1, phi1, acr, fcr] = estimator1(frame, fs, w, nfft, th)
    # freq1 = omega1 / (2 * np.pi * fs)
    # parameter = Parameters(a1, freq1, phi1)
    parameter = params.get_parameter(i)
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()
    vector_frames[i, :] = s

# Define an empty vector to receive result
vector_time = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)
s_gt_env = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)

# Loop for each frame and overlap-add
time_index = 0
for i in xrange(number_frames):
    vector_time[time_index:time_index + window_size] += vector_frames[i, :]
    s_gt_env[time_index:time_index + window_size] += w
    time_index += synthesis_hop

# Comparison
plt.title("Synthetised and original frame")
plt.plot(range(len(vector_time)), vector_time, range(len(vector_time)), s_gt_env*s_gt[:len(vector_time)], '--')
plt.xlabel("Echantillons")
plt.ylabel("Signal")
print "RMSE = " + str((1./len(vector_time))*np.sqrt(np.sum((vector_time-s_gt[:len(vector_time)])**2)))

plt.show()

