#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from SpectrumGenerator import *
from NonStationaryLUT import *  # todo : Non stationary LUT


class NonStationarySpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class NonStationarySpectrumGenerator(SpectrumGenerator):
    def __init__(self,regular_lut,_window_size, parameters,  acr_domain, fcr_domain):
        self._regular_lut = regular_lut
        self._lut = NonStationaryLUT(regular_lut, acr_domain, fcr_domain)
        SpectrumGenerator.__init__(self,_window_size, parameters)

    def _add_lobe(self, k):  # todo : Implement _add_lobe for non stationary signals
        pass
