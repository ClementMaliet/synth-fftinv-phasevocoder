#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET


# Import the libraries

from scipy import signal
from data_structure import *


class StationaryLobe(object):
    def __init__(self, window_type, window_size, nfft):
        self._window_type = window_type
        self._window_size = window_size
        self._nfft = nfft
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._lobe = Spectrum([], [])
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

        w1 = np.fft.fft(self._window, self._window_size + self._nfft)

        # Correct phase
        #sw = np.zeros(self._nfft)
        #sw[:(self._window_size - 1) / 2] = self._window[(self._window_size - 1) / 2:]
        #sw[self._nfft - (self._window_size - 1) / 2:] = s[:(self._window_size + 1) / 2]

        # Compute the fft
        #w1 = np.fft.fft(sw, self._nfft)
        #mod_fft_s = np.absolute(w1)

        w1 = np.concatenate([w1[-4:], w1[:5]])
        self._lobe.set_complex_spectrum(w1)

    def get_lobe(self):
        return self._lobe
