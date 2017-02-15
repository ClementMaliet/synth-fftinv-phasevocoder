from abc import ABCMeta, abstractmethod

from data_structure import *


class SpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass


class SpectrumGenerator:

    __metaclass__ = ABCMeta

    """The class SpectrumGenerator is used to generate a spectrum with 9 points per lobe"""

    def __init__(self, _window_size, parameters, spectrum):
        self._parameters = parameters
        self._spectrum = spectrum
        self._window_size = _window_size

    @abstractmethod
    def _add_lobe(self, k):
        pass

    def set_parameters(self, new_parameters):
        self._spectrum = Spectrum.void_spectrum(self._window_size)
        self._parameters = new_parameters

    def get_spectrum(self):
        for k in range(self._parameters.get_number_sinuses()):
            self._add_lobe(k)