#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from data_structure import *
from abc import ABCMeta, abstractmethod
from scipy import signal


class LobeGenerator(object):
    __metaclass__ = ABCMeta

    def __init__(self, window_type, window_size, nfft):
        self._window_type = window_type
        self._window_size = window_size
        self._nfft = nfft
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= np.sum(self._window)
        self._lobe = Spectrum([], [])

    def _set_window_size(self, window_size):
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lobe()
    window_size = property(fset=_set_window_size)

    def _set_window_type(self, window_type):
        self._window_type = window_type
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lobe()
    window_type = property(fset=_set_window_type)

    @abstractmethod
    def _gen_lobe(self):
        pass

    def get_lobe(self):
        return self._lobe
