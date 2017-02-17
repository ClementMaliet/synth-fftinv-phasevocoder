from sources.Synthesizer import *
from scipy import signal
import matplotlib.pyplot as plt

parameter = Parameters(np.array([1]), np.array([0.12]), np.array([1]))
print 1
synth = StationarySynthesizer(1015, "hamming", 3, 500, 500, parameter)
print 2
synth.set_next_frame(parameter)
print 3
s = synth.get_next_frame()
print 4
w = signal.get_window("hamming", 1015)
s_gt = np.array([1])*np.exp(2*np.pi*1j*0.12*np.array(range(1015)) + 1j)
s_gt *= w
print 5
plt.plot(range(1015), s_gt, range(1015), s)
plt.xlabel("Echantillons")
plt.ylabel("Signal")
plt.legend("Signal de reference", "Signal synthetise")
plt.show()