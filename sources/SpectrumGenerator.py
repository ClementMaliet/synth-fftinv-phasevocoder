#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from abc import ABCMeta, abstractmethod

from data_structure import *


class SpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class SpectrumGenerator:

    __metaclass__ = ABCMeta

    """The class SpectrumGenerator is used to generate a spectrum with 9 points per lobe"""

    def __init__(self, window_size, parameters):
        self._parameters = parameters
        self._spectrum = Spectrum.void_spectrum(window_size)
        self._window_size = window_size

    @abstractmethod
    def _add_lobe(self, k):
        pass

    def _set_parameters(self, new_parameters):
        self._spectrum = Spectrum.void_spectrum(self._window_size)
        self._parameters = new_parameters
    parameters = property(fset=_set_parameters)

    def get_spectrum(self):
        for k in range(self._parameters.get_number_sinuses()):
            self._add_lobe(k)
        return self._spectrum
