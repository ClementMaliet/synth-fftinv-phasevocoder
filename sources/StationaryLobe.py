#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET


# Import the libraries

from scipy import signal
import numpy as np
from data_structure import *



# Classes:

class StationaryLobe:
    def __init__(self, window_type, window_size):
        self._window_type = window_type
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._lobe = Spectrum([],[])
        self._gen_lobe()

    def _set_window_size(self, window_size):
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lobe()

    def _set_window_type(self, window_type):
        self._window_type = window_type
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lobe()
    window_type = property(fset=_set_window_type)
    window_size = property(fset=_set_window_size)

    def _gen_lobe(self):

        w1 = np.fft.fft(self._window, self._window_size)
        w1 = np.concatenate([w1[-4:], w1[:5]])
        self._lobe.set_complex_spectrum(w1)



    def get_lobe(self):
        return self._lobe
