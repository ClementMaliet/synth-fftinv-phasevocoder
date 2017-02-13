#!/usr/bin/python
# -*- coding: utf-8 -*-

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
        if np.iscomplexobj(complex_spectrum) is False:
            raise InconsistentSpectrumError('Real spectrum cannot be converted to amplitude and phase')
        else:
            cls._amplitude = np.absolute(complex_spectrum)
            cls._phase = np.angle(complex_spectrum)
            cls._nfft = len(complex_spectrum)

    @classmethod
    def void_spectrum(cls, nfft):
        cls._amplitude = np.zeros(nfft)
        cls._phase = np.zeros(nfft)
        cls._nfft = nfft

    def set_spectrum(self, amplitude, phase):
        if len(amplitude) != len(phase):
            raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
        else:
            self._amplitude = amplitude
            self._phase = phase
            self._nfft = len(amplitude)

    def set_complex_spectrum(self, complex_spectrum):
        if np.iscomplexobj(complex_spectrum) is False:
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

    def get_nfft(self):
        return self._nfft


class ParametersError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class InconsistentParametersError(ParametersError):
    """Base class for inconsistency of parameters, raised when the numbers of sinusoid and/or parameters do not
    add up."""

    def __init__(self, message):
        self.message = message


class NegativeAmplitudeError(InconsistentParametersError):
    """Raised if one of the amplitude provided is negative"""


class BoundParametersError(ParametersError):
    """Raised if asked for a non existing sinusoid's parameters (index out of bound)
    Attributes:
            message -- explanation of the error
        """

    def __init__(self, message):
        self.message = message


class Parameters:
    """The class parameters is used to store the parameters of the model's stationary sinusoid"""

    def __init__(self, amplitudes, frequencies, phases):
        if len(amplitudes) != len(frequencies) or len(amplitudes) != len(phases) or len(frequencies) != len(phases):
            raise InconsistentParametersError("Sinusoid parameters provided are not the same length")
        elif any([a < 0 for a in amplitudes]):
            raise NegativeAmplitudeError("One of the amplitude provided is negative")
        else:
            self._amplitudes = amplitudes
            self._frequencies = frequencies
            self._phases = phases
            self._number_sinuses = len(amplitudes)

    def get_amplitude(self, k):
        if k < 0:
            raise BoundParametersError("Index negative !")
        elif k > self._number_sinuses - 1:
            raise BoundParametersError("Index out of parameters range")
        else:
            return self._amplitudes[k]

    def get_frequency(self, k):
        if k < 0:
            raise BoundParametersError("Index negative !")
        elif k > self._number_sinuses - 1:
            raise BoundParametersError("Index out of parameters range")
        else:
            return self._frequencies[k]

    def get_phase(self, k):
        if k < 0:
            raise BoundParametersError("Index negative !")
        elif k > self._number_sinuses - 1:
            raise BoundParametersError("Index out of parameters range")
        else:
            return self._phases[k]


class NonStationaryParameters(Parameters):
    """Subclass of Parameters which yields 2 extra parameters and their accessors :
        - acr : Amplitude Change Rate
        - fcr : Frequency Change Rate"""

    def __init__(self, amplitudes, frequencies, phases, acrs, fcrs):
        Parameters.__init__(amplitudes, frequencies, phases)
        if len(acrs) != len(fcrs) or len(acrs) != len(amplitudes) or len(fcrs) != len(amplitudes):
            raise InconsistentParametersError("Sinusoid parameters provided are not the same length")
        else:
            self._acrs = acrs
            self._fcrs = fcrs

    def get_acrs(self, k):
        if k < 0:
            raise BoundParametersError("Index negative !")
        elif k > self._number_sinuses - 1:
            raise BoundParametersError("Index out of parameters range")
        else:
            return self._acrs[k]

    def get_fcrs(self, k):
        if k < 0:
            raise BoundParametersError("Index negative !")
        elif k > self._number_sinuses - 1:
            raise BoundParametersError("Index out of parameters range")
        else:
            return self._fcrs[k]