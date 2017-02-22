from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
sine_duration = 1
window_length = 23e-3

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hamming"
zero_padding_factor = 2
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 500
synthesis_hop = 500

parameter = Parameters(np.array([1]), np.array([0.1]), np.array([3.]))

# Find the max number of slices that can be obtained
number_frames = np.floor((sine_duration*fs-window_size)/analysis_hop)

# Create a matrix with time slices
vector_frames = np.zeros((number_frames, window_size))

# Analysis
w = signal.get_window(window_type, window_size)
w /= np.sum(w)

s_gt = np.zeros(sine_duration*fs)
for i in xrange(parameter.get_number_sinuses()):
    s_gt += parameter.get_amplitude(i)*np.exp(2*np.pi*1j*parameter.get_frequency(i)*np.array(range(sine_duration*fs)) +
                                              parameter.get_phase(i)*1j).real

for i in xrange(number_frames):
    vector_frames[i, :] = s_gt[i*window_size:(i+1)*window_size]
    vector_frames[i, :] *= w

    sw_gt = np.zeros(nfft)
    sw_gt[:(window_size - 1) / 2] = vector_frames[i, ((window_size + 1) / 2):]
    sw_gt[nfft - (window_size - 1) / 2 - 1:] = vector_frames[i, :(window_size + 1) / 2]

    sfft = np.fft.fft(sw_gt, nfft)
    phases = np.array([np.angle(sfft[int(round(a * nfft))]) for a in parameter._frequencies])

    # Synthesis
    parameter._phases = phases
    synth = StationarySynthesizer(window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop, parameter)
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()

# Define an empty vector to receive result
vector_time = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size, 1);

# Loop for each frame and overlap-add todo : UNFINISHED overlap add, look at matlab code below in commentary
time_index = 1
for i in xrange(number_frames)

    vector_time[time_index:time_index + window_size - 1] =

time_index = time_index + synthesis_hop



# Comparison
omega = np.arange(nfft) / float(nfft) - 0.5
sfft = np.fft.fftshift(sfft)

plt.figure()
plt.title("Synthetized spectrum")
plt.subplot(2,1,1)
plt.plot(omega, 10*np.log10(np.fft.fftshift(synth._current_spectrum._amplitude)**2), omega, 10*np.log10(np.absolute(sfft)**2))
plt.subplot(2,1,2)
plt.plot(omega, np.fft.fftshift(synth._current_spectrum._phase), omega, 2*np.angle(sfft))
plt.figure()
plt.title("Synthetised and original frame")
plt.plot(range(window_size), s, range(window_size), s_gt)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
print "RMSE = " + str((1./window_size)*np.sqrt(np.sum((s-s_gt)**2)))

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