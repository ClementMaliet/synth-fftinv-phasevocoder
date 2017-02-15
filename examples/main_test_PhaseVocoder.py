from sources.SpectrumGenerator import SpectrumGenerator
from sources.StationaryLobe import *
from sources.PhaseVocoder import *
from scipy import signal
import matplotlib.pyplot as plt

from sources.StationarySpectrumGenerator import StationarySpectrumGenerator

window_size = 1024
window_type = 'hamming'

s = StationaryLobe(window_type, window_size)

parameters = Parameters(np.array([2, 3]), np.array([0.3, 0.4]), np.array([50, 40]))

sg = StationarySpectrumGenerator(window_type, window_size, parameters)

pv_stat_test = StationaryPhaseVocoder(1, 4, sg.get_spectrum()) #(Analyse_hop, Synthesis_hop, Current_analysis_spectrum)

#Synthetic Spectrum
plt.figure()
plt.plot(np.linspace(0, 1, window_size), sg.get_spectrum()._amplitude)

#Phase Vocoded Spectrum
plt.figure()
plt.plot(np.linspace(0, 1, window_size), pv_stat_test.get_pv_spectrum()._amplitude)
plt.show()
