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
    def from_complex_spectrum(cls, complex_spectrum, x):
        if np.iscomplexobj(complex_spectrum) is False:
            raise InconsistentSpectrumError('Real spectrum cannot be converted to amplitude and phase')
        else:
            amplitude = np.absolute(complex_spectrum)
            phase = np.angle(complex_spectrum)
            return cls(amplitude, phase, x)

    @classmethod
    def void_spectrum(cls, nfft):
        amplitude = np.zeros(nfft)
        phase = np.zeros(nfft)
        x = np.zeros(nfft)
        return cls(amplitude, phase, x)

    def set_spectrum(self, amplitude, phase, start_bin=None, stop_bin=None, x):
        Spectrum.set_spectrum(self, amplitude, phase, start_bin=None, stop_bin=None)
        self._x = x

    def set_complex_spectrum(self, complex_spectrum, start_bin=None, stop_bin=None, x):
        Spectrum.set_complex_spectrum(self,complex_spectrum,start_bin=None, stop_bine=None)
        self._x = x

    def get_amplitude(self, k):
        if k < 0:
            raise BoundSpectrumError("Index negative !")
        elif k > self._nfft - 1:
            raise BoundSpectrumError("Index out of spectrum range")
        else:
            return self._amplitude[k]

    def get_phase(self, k):
        if k < 0:
            raise BoundSpectrumError("Index negative !")
        elif k > self._nfft - 1:
            raise BoundSpectrumError("Index out of spectrum range")
        else:
            return self._phase[k]

    def get_abscissa(self, k):
        return self._x[k]