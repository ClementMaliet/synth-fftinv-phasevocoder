from sources.SpectrumGenerator import SpectrumGenerator
from sources.StationaryLobeGenerator import *
from sources.PhaseVocoder import *
from scipy import signal
import matplotlib.pyplot as plt

from sources.StationarySpectrumGenerator import StationarySpectrumGenerator

window_size = 1025
window_type = 'hamming'
nfft = 8192

s = StationaryLobeGenerator(window_type, window_size, nfft)

parameters = Parameters(np.array([2, 3]), np.array([0.3, 0.4]), np.array([50, 40]))

sg = StationarySpectrumGenerator(window_type, window_size, parameters, nfft)

pv_stat_test = StationaryPhaseVocoder(1, 4, sg.get_spectrum()) #(Analyse_hop, Synthesis_hop, Current_analysis_spectrum)

#Synthetic Spectrum
plt.figure()
plt.plot(np.linspace(0, 1, nfft), sg.get_spectrum()._amplitude)
plt.title("synthesised spectrum")
#Phase Vocoded Spectrum
plt.figure()
plt.plot(np.linspace(0, 1, nfft), pv_stat_test.get_pv_spectrum()._amplitude)
plt.title("Phase vocoded spectrum")
plt.show()

plt.close()