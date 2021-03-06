#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

# This file contains the definition of the Spectrum and the Parameters class with their respective subclasses
# The Spectrum class is the class used to exchange and store full spectra and main lobes. The spectrum points are stored
# as amplitude and unwrapped phase and are accessed in a point wise fashion.
# The Parameters class is the class used to exchange and store the sinusoid model parameters which consist of :
#     - alpha, f and phi for stationary sinusoid
#     - alpha, mu, f, psi and phi for non stationary sinusoid
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


class AddSpectrumError(SpectrumError):
    """Raised if asked to add a spectrum with a compatible but non viable element
    Attributes:
            message -- explanation of the error
        """

    def __init__(self, message):
        self.message = message


class Spectrum(object):
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
            amplitude = np.absolute(complex_spectrum)
            phase = np.angle(complex_spectrum)
            return cls(amplitude, phase)

    @classmethod
    def void_spectrum(cls, nfft):
        amplitude = np.zeros(nfft)
        phase = np.zeros(nfft)
        return cls(amplitude, phase)

    def __add__(self, other):
        if isinstance(other, self.__class__):
            if other._nfft == self._nfft:
                return Spectrum.from_complex_spectrum(self._amplitude*np.exp(1j*self._phase) +
                                                      other._amplitude*np.exp(1j*other._phase))
            else:
                raise AddSpectrumError("Addition of non equally sized spectrum attempted")
        elif isinstance(other, np.ndarray) and np.iscomplexobj(other):
            if len(other) == self._nfft:
                return Spectrum.from_complex_spectrum(self._amplitude*np.exp(1j*self._phase) + other)
            else:
                raise AddSpectrumError("Addition of non equally sized spectrum attempted")
        else:
                raise NotImplementedError
    __radd__ = __add__

    def __iadd__(self, other):
        if isinstance(other, self.__class__):
            if other._nfft == self._nfft:
                self.set_complex_spectrum(self._amplitude*np.exp(1j*self._phase) +
                                          other._amplitude*np.exp(1j*other._phase))
            else:
                raise AddSpectrumError("Addition of non equally sized spectrum attempted")
        elif isinstance(other, np.ndarray) and np.iscomplexobj(other):
            if len(other) == self._nfft:
                self.set_complex_spectrum(self._amplitude * np.exp(1j * self._phase) + other)
            else:
                raise AddSpectrumError("Addition of non equally sized spectrum attempted")
        else:
            raise NotImplementedError
        return self

    def __mul__(self, other):  # todo : multiplication with a single non array scalar
        if isinstance(other, np.ndarray) :
            return Spectrum.from_complex_spectrum(self._amplitude*np.exp(1j*self._phase) * other)
        else:
            raise NotImplementedError
    __rmul__ = __mul__

    def __imul__(self, other):  # todo : inline multiplication with a single non array scalar
        if isinstance(other, np.ndarray):
            self.set_complex_spectrum(self._amplitude*np.exp(1j*self._phase) * other)
        else:
            raise NotImplementedError
        return self

    def set_spectrum(self, amplitude, phase, start_bin=None, stop_bin=None):
        if start_bin is None and stop_bin is None:
            if len(amplitude) != len(phase):
                raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
            else:
                self._amplitude = amplitude
                self._phase = phase
                self._nfft = len(amplitude)
        elif start_bin is None and stop_bin is not None:
            if len(amplitude) != len(phase):
                raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
            elif stop_bin > self._nfft:
                raise BoundSpectrumError("Attempted to modify non existing bins (stop_bin too large)")
            else:
                self._amplitude[:stop_bin] = amplitude
                self._phase[:stop_bin] = phase
        elif start_bin is not None and stop_bin is None:
            if len(amplitude) != len(phase):
                raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
            elif start_bin < 0:
                raise BoundSpectrumError("Attempted to modify non existing bins (start_bin negative)")
            else:
                self._amplitude[:stop_bin] = amplitude
                self._phase[:stop_bin] = phase
        elif start_bin is not None and stop_bin is not None:
            if len(amplitude) != len(phase):
                raise InconsistentSpectrumError('Amplitude and Phase provided have different length')
            elif stop_bin > self._nfft or start_bin < 0:
                raise BoundSpectrumError("Attempted to modify non existing bins (interval incorrect)")
            else:
                self._amplitude[start_bin:stop_bin] = amplitude
                self._phase[start_bin:stop_bin] = phase

    def set_complex_spectrum(self, complex_spectrum, start_bin=None, stop_bin=None):
        if start_bin is None and stop_bin is None:
            if not np.iscomplexobj(complex_spectrum):
                raise InconsistentSpectrumError('Spectrum is not complex')
            else:
                self._amplitude = np.absolute(complex_spectrum)
                self._phase = np.angle(complex_spectrum)
                self._nfft = len(complex_spectrum)
        elif start_bin is None and stop_bin is not None:
            if not np.iscomplexobj(complex_spectrum):
                raise InconsistentSpectrumError('Spectrum is not complex')
            elif stop_bin > self._nfft:
                raise BoundSpectrumError("Attempted to modify non existing bins (stop_bin too large)")
            else:
                self._amplitude[:stop_bin] = np.absolute(complex_spectrum)
                self._phase[:stop_bin] = np.angle(complex_spectrum)
        elif start_bin is not None and stop_bin is None:
            if not np.iscomplexobj(complex_spectrum):
                raise InconsistentSpectrumError('Spectrum is not complex')
            elif start_bin < 0:
                raise BoundSpectrumError("Attempted to modify non existing bins (start_bin negative)")
            else:
                self._amplitude[:stop_bin] = np.absolute(complex_spectrum)
                self._phase[:stop_bin] = np.angle(complex_spectrum)
        elif start_bin is not None and stop_bin is not None:
            if not np.iscomplexobj(complex_spectrum):
                raise InconsistentSpectrumError('Spectrum is not complex')
            elif stop_bin > self._nfft or start_bin < 0:
                raise BoundSpectrumError("Attempted to modify non existing bins (interval incorrect)")
            else:
                self._amplitude[start_bin:stop_bin] = np.absolute(complex_spectrum)
                self._phase[start_bin:stop_bin] = np.angle(complex_spectrum)

    def get_complex_spectrum(self):
        return self._amplitude*np.exp(1j*self._phase)

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

    def _get_full_amplitude(self):
        return self._amplitude
    amplitude = property(fget=_get_full_amplitude)

    def _get_full_phase(self):
        return self._phase
    phase = property(fget=_get_full_phase)


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


class Parameters(object):
    """The class parameters is used to store the parameters of the model's stationary sinusoid"""

    def __init__(self, amplitudes, frequencies, phases, regions=None, peak_indices=None):
        if len(amplitudes) != len(frequencies) or len(amplitudes) != len(phases) or len(frequencies) != len(phases):
            raise InconsistentParametersError("Sinusoid parameters provided are not the same length")
        elif any([a < 0 for a in amplitudes]):
            raise NegativeAmplitudeError("One of the amplitude provided is negative")
        else:
            self._amplitudes = amplitudes
            self._frequencies = frequencies
            self._phases = phases
            self._number_sinuses = len(amplitudes)
            self._regions = regions if regions is not None else np.array([0, self._number_sinuses])
            self._peak_indices = peak_indices if peak_indices is not None else np.arange(self._number_sinuses)

    @classmethod
    def void_parameters(cls, number_sinuses):
        amplitudes = np.zeros(number_sinuses)
        phases = np.zeros(number_sinuses)
        frequencies = np.zeros(number_sinuses)
        return cls(amplitudes, frequencies, phases)

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

    def get_number_sinuses(self):
        return self._number_sinuses

    def get_regions(self):
        # set the min amplitude of a peak
        min_amp = 1
        if self._amplitudes.max > min_amp:
            self._peak_indices = self._amplitudes.find_peaks_cwt()

            # to attache a sample to a peak we detect the minimum between each peak
            self._regions[0] = 0
            previous_peak = self._peak_indices(0)
            for k in range(len(self._peak_indices) - 1):
                self._regions[k+1] = np.amin(self._amplitudes[previous_peak:self._peak_indices(k+1)])
                previous_peak =  self._peak_indices(k+1)

            self._regions[len( self._peak_indices)+1] = len(self._amplitudes)
        # if max of amplitude is too short, we consider the whole frame as a regions of interest

        else:
            self._peak_indices = self._amplitudes.max
            self._regions[0] = 0
            self._regions[1] = len(self._amplitudes)
        return self._peak_indices, self._regions


class NonStationaryParameters(Parameters):
    """Subclass of Parameters which yields 2 extra parameters and their accessors :
        - acr : Amplitude Change Rate
        - fcr : Frequency Change Rate"""

    def __init__(self, amplitudes, frequencies, phases, acrs, fcrs):
        Parameters.__init__(self, amplitudes, frequencies, phases)
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
