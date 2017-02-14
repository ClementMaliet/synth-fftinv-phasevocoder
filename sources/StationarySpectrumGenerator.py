#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

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



