from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

fs = 44100
sine_duration = 46e-3
window_length = 23e-3

window_size = int(round(fs*window_length) if round(fs*window_length) % 2 != 0 else round(fs*window_length) + 1)
window_type = "hanning"
zero_padding_factor = 3
nfft = 2**(next2pow(window_size) + zero_padding_factor)
analysis_hop = 505
synthesis_hop = 505

parameter = Parameters(np.array([1]), np.array([0.1]), np.array([3.]))

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
                                              np.array(range(int(round(sine_duration*fs)))) +
                                              parameter.get_phase(i)*1j).real

for i in xrange(number_frames):
    # Synthesis
    synth.set_next_frame(parameter)
    s = synth.get_next_frame()
    vector_frames[i, :] = s

# Define an empty vector to receive result
vector_time = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)
s_gt_env = np.zeros(number_frames * synthesis_hop - synthesis_hop + window_size)

# Loop for each frame and overlap-add todo : UNFINISHED overlap add, look at matlab code below in commentary
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