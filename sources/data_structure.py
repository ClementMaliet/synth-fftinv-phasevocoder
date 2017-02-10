#!/usr/bin/python
# -*- coding: utf8-*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

# This file contains the definition of the Spectrum and the Parameters class with their respective subclasses
# The Spectrum class is the class used to exchange and store full spectra and main lobes. The spectrum points are stored
# as amplitude and unwrapped phase and are accessed in a point wise fashion.
# The Parameters class is the class used to exchange and store the sinusoid model parameters which consist of :
#     - lambda, omega and phi for stationary sinusoid
#     - lambda, mu, omega, psi, phi for non stationary sinusoid
# They are then accessed in a sinus wise fashion.

import numpy as np


class SpectrumError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class InconsistentSpectrumError(SpectrumError):
    """Raised if the spectrum provided is inconsistent and cannot be a valid spectrum are converted to a valid spectrum
    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class BoundSpectrumError(SpectrumError):
    """Raised if asked for a non existing point of the spectrum (index out of bound)
    Attributes:
            message -- explanation of the error
        """

    def __init__(self, message):
        self.message = message


class Spectrum:
    """The class Spectrum is used to store a spectrum with amplitude and phase."""
    def __init__(self, amplitude, phase):
        if len(amplitude) != len(phase):
            raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
        else:
            self._amplitude = amplitude
            self._phase = phase
            self._nfft = len(amplitude)

    @classmethod
    def from_complex_spectrum(cls, complex_spectrum):
        if type(complex_spectrum) != "complex_":
            raise InconsistentSpectrumError('Real spectrum cannot be converted to amplitude and phase')
        else:
            cls._amplitude = np.absolute(complex_spectrum)
            cls._phase = np.angle(complex_spectrum)
            cls._nfft = len(complex_spectrum)

    def set_spectrum(self, amplitude, phase):
        if len(amplitude) != len(phase):
            raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
        else:
            self._amplitude = amplitude
            self._phase = phase
            self._nfft = len(amplitude)

    def set_complex_spectrum(self, complex_spectrum):
        if type(complex_spectrum) != np.complex_:
            raise InconsistentSpectrumError('Real spectrum cannot be converted to amplitude and phase')
        else:
            self._amplitude = np.absolute(complex_spectrum)
            self._phase = np.angle(complex_spectrum)
            self._nfft = len(complex_spectrum)

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
