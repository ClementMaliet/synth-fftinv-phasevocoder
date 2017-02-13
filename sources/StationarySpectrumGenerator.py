# Import :

from SpectrumGenerator import *
from StationaryLobe import *

# Classes

class StationarySpectrumGenerator(SpectrumGenerator) :

    def __init__(self, window_type,window_size, parameters):
        SpectrumGenerator.init(self, window_size, parameters)
        self._stationary_lobe = StationaryLobe(window_type, window_size)

    def _add_lobe(self, k):
        lobe = StationaryLobe.get_lobe()
        amplitude = self.parameters.get_amplitude(k)
        phase = self.parameters.get_phase(k)
        frequency = self.parameters.get_frequency(k)



