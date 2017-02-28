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

# Find the max number of slices that can be obtained
number_frames = int(np.floor((sine_duration*fs-window_size)/analysis_hop))

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)


for i in xrange(number_frames):
    vector_frames[i, :] = s_gt[i*analysis_hop:i*analysis_hop + window_size]
    #vector_frames[i, :] *= w

    [a1, omega1, phi1, acr, fcr] = estimator1(vector_frames[i, :], fs, w, nfft, th)
    freq1 = omega1 /(2*np.pi*fs)
    parameter = Parameters(a1, freq1, phi1)

    #parameter = Parameters(np.array([20939]), np.array([3671.09384348]), np.array([-1.7773237]))

    sw_gt = np.zeros(nfft)
    sw_gt[:(window_size - 1) / 2] = vector_frames[i, ((window_size + 1) / 2):]
    sw_gt[nfft - (window_size - 1) / 2 - 1:] = vector_frames[i, :(window_size + 1) / 2]

    sfft = np.fft.fft(sw_gt, nfft)
    # phases = np.array([np.angle(sfft[int(round(a * nfft))]) for a in parameter._frequencies])

    # Synthesis
    synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)

    # parameter._phases = phases
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

