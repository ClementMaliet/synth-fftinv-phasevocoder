#!/usr/bin/python
# -*- coding: utf-8 -*-

# Non stationary lobe with amplitude, phase and frequencies
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET


#Import the libraries:
from data_structure import *


#Class:

class NonStationaryLobe(Spectrum):

    def __init__(self, amplitude, phase, x):
        Spectrum.__init__(self, amplitude, phase)
        self._x = x

    @classmethod
    def from_complex_lobe(cls, complex_spectrum, x):
        if np.iscomplexobj(complex_spectrum) is False:
            raise InconsistentSpectrumError('Real spectrum cannot be converted to amplitude and phase')
        else:
            amplitude = np.absolute(complex_spectrum)
            phase = np.angle(complex_spectrum)
            return cls(amplitude, phase, x)

    @classmethod
    def void_lobe(cls, nfft):
        amplitude = np.zeros(nfft)
        phase = np.zeros(nfft)
        x = np.zeros(nfft)
        return cls(amplitude, phase, x)

    def set_lobe(self, amplitude, phase, x, start_bin=None, stop_bin=None):
        Spectrum.set_spectrum(self, amplitude, phase, start_bin=None, stop_bin=None)
        self._x = x

    def set_complex_lobe(self, complex_spectrum, x, start_bin=None, stop_bin=None):
        Spectrum.set_complex_spectrum(self, complex_spectrum, start_bin=None, stop_bin=None)
        self._x = x

    def get_abscissa(self, k):

        return self._x[k]