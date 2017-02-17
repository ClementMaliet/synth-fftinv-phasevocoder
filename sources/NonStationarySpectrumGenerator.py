#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from SpectrumGenerator import *
from NonStationaryLUT import *  # todo : LobeGeneration in NonStationaryLUT


class NonStationarySpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class NonStationarySpectrumGenerator(SpectrumGenerator):
    def __init__(self, window_type, window_size, parameters, regular_grid, acr_domain,
                 fcr_domain, number_acr, number_fcr):
        self._regular_lut = regular_grid
        self._lobe_generator = NonStationaryLUT(regular_grid, acr_domain, fcr_domain, number_acr,
                                                number_fcr, window_type, window_size)
        SpectrumGenerator.__init__(self, window_type, window_size, parameters)

    def _add_lobe(self, k):  # todo : Implement _add_lobe for non stationary signals
        pass
