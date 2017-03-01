from sources.Synthesizer import *
import matplotlib.pyplot as plt
from estimation_func import *
from numpy import *
from matplotlib.pyplot import *
from scipy import signal
from scipy.io import wavfile

filename = "Mix.wav"
[fs, x] = wavfile.read(filename)

sine_duration = len(x)/fs  # s
window_length = 23e-3  # s

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 0
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 505
synthesis_hop = 505
th = -10

s_gt = x

synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                              Parameters.void_parameters(1))
# Find the max number of slices that can be obtained
number_frames = int(np.floor((sine_duration*fs-window_size)/analysis_hop))

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

# We compute the true initial phase to take into account the windowing effect on phase
frame = np.zeros(window_size)
frame[:window_size] = s_gt[:window_size]
frame *= w

sw_gt = np.zeros(nfft)
sw_gt[:(window_size - 1) / 2] = frame[((window_size + 1) / 2):]
sw_gt[nfft - (window_size - 1) / 2 - 1:] = frame[:(window_size + 1) / 2]

[a1, omega1, phi1, acr, fcr] = estimator1(frame, fs, w, nfft, th)
freq1 = omega1 /(2*np.pi*fs)
parameter = Parameters(a1, freq1, phi1)

for i in xrange(number_frames):
    frame[:window_size] = s_gt[i*analysis_hop:i*analysis_hop + window_size]

    [a1, omega1, phi1, acr, fcr] = estimator1(frame, fs, w, nfft, th)
    freq1 = omega1 / (2 * np.pi * fs)
    parameter = Parameters(a1, freq1, phi1)

    parameter._phases = np.array([synth.current_spectrum.phase[int(round(a * nfft))] for a in parameter._frequencies])
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

