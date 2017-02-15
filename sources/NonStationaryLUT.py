#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

import numpy as np
from scipy import signal
import warnings


class NonStationaryLUT(object):
    def __init__(self, regular_grid, acr_domain, fcr_domain, window_type, window_size):
        self._regular_grid = regular_grid
        self._domain = [acr_domain, fcr_domain]
        self._window_type = window_type
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= np.sum(self._window)  # We normalize the window
        self._gen_lut()

    def _gen_lut(self):
        if self._regular_grid:
            self._gen_uniform_lut()
        else:
            self._gen_non_uniform_lut()

    def _set_window_size(self, window_size):
        warnings.warn("LUT will be recomputed, it may take a while", UserWarning)
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lut()

    def _set_window_type(self, window_type):
        warnings.warn("LUT will be recomputed, it may take a while", UserWarning)
        self._window_type = window_type
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self._gen_lut()
    window_type = property(fset=_set_window_type)
    window_size = property(fset=_set_window_size)

    def _gen_uniform_lut(self):
        pass

    def _gen_non_uniform_lut(self):
        pass
