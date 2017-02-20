#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from abc import ABCMeta, abstractmethod

from data_structure import *
from LobeGenerator import LobeGenerator


class SpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class SpectrumGenerator(object):

    __metaclass__ = ABCMeta

    """The class SpectrumGenerator is used to generate a spectrum with 9*self._zero_padding_factor points per lobe"""

    def __init__(self, window_size, parameters, nfft):
        self._parameters = parameters
        self._nfft = nfft
        self._spectrum = Spectrum.void_spectrum(self._nfft)
        self._window_size = window_size
        # The structure imposes :
        # self._lobe_generator = LobeGenerator(window_type, window_size, nfft)

    @abstractmethod
    def _add_lobe(self, k):
        pass

    def _set_window_size(self, window_size):
        self._window_size = window_size
        self._lobe_generator.window_size = window_size
    window_size = property(fset=_set_window_size)

    def _set_window_type(self, window_type):
        self._window_type = window_type
        self._lobe_generator.window_type = window_type
    window_type = property(fset=_set_window_type)

    def _set_parameters(self, new_parameters):
        print "Set parameters"
        self._spectrum = Spectrum.void_spectrum(self._nfft)
        self._parameters = new_parameters
    parameters = property(fset=_set_parameters)

    def get_spectrum(self):
        print "Make spectrum"
        for k in range(self._parameters.get_number_sinuses()):
            self._add_lobe(k)
        return self._spectrum
