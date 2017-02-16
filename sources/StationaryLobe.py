#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET


# Import the libraries

from LobeGenerator import LobeGenerator
import numpy as np


class StationaryLobe(LobeGenerator):
    def __init__(self, window_type, window_size, nfft):
        LobeGenerator.__init__(self, window_type, window_size, nfft)
        self._gen_lobe()

    def _gen_lobe(self):
        w1 = np.fft.fft(self._window, self._window_size)
        w1 = np.concatenate([w1[-4:], w1[:5]])
        self._lobe.set_complex_spectrum(w1)
