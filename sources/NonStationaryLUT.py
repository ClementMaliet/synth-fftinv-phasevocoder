#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from scipy import signal
from scipy.interpolate import interp1d
import warnings
from data_structure import *
from Synthesizer import next2pow


class NonStationaryLUT(LobeGenerator):  # todo : implement LobeGenerator to use polymorphism
    def __init__(self, regular_grid, acr_domain, fcr_domain, number_acr, number_fcr, window_type, window_size, fs=None):
        if fs is None:
            self._fs = 44100
        else:
            self._fs = fs
        self._regular_grid = regular_grid
        self._domain = [acr_domain, fcr_domain]
        self._number_acr = number_acr
        self._number_fcr = number_fcr
        self._number_points = number_acr*number_fcr
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
        acr = np.linspace(self._domain[0][0], self._domain[0][1], self._number_acr)
        fcr = np.linspace(self._domain[1][0], self._domain[1][1], self._number_fcr)
        self._sample_grid = [(a, b) for a in acr for b in fcr]
        n = np.array(range(self._window_size))
        t = (1/self._fs)*n
        zero_padding_factor = 9
        nfft = 2**(next2pow(self._window_size)+zero_padding_factor)

        for i in xrange(self._number_acr):
            for j in xrange(self._number_fcr):
                pass

    def _gen_non_uniform_lut(self):
        pass

    def _gen_lobe_legacy(self, i, j, acr, fcr, t, n, nfft):
        # Generate chirp
        mu = acr[i]
        psi = fcr[j]
        s = (1 + mu * t) * np.exp(1j * ((0.5 * psi * t ** 2) + (2 * np.pi * 0.0227 * n)))

        # Apply window
        s *= self._window

        # Correct phase
        sw = np.zeros(nfft)
        sw[:(self._window_size - 1) / 2] = s[(self._window_size - 1) / 2:]
        sw[nfft - (self._window_size - 1) / 2:] = s[:(self._window_size + 1) / 2]

        # Compute the fft
        fft_s = np.fft.fft(sw, nfft)
        mod_fft_s = np.absolute(fft_s)

        # Search for the main lobe
        index = np.argmax(mod_fft_s)

        # Initialize tests to find the zeros
        left_peak = mod_fft_s[:-2] > mod_fft_s[1:-1]
        right_peak = mod_fft_s[2:] > mod_fft_s[1:-1]

        # Find location of zeros in the FFT
        mod_zeros_loc = np.nonzero(np.logical_and(left_peak, right_peak))[0]

        # Find the closest zero to the lob peak
        index_min = np.argmin(np.absolute(index - mod_zeros_loc))

        # Set the lower and upper zero of the lobe
        if mod_zeros_loc[index_min] > index:
            lower_zero_loc = mod_zeros_loc[index_min - 1]
            upper_zero_loc = mod_zeros_loc[index_min]
        else:
            lower_zero_loc = mod_zeros_loc[index_min]
            upper_zero_loc = mod_zeros_loc[index_min + 1]

        # Indexes of the main lobe
        lobe_index = np.array(range(lower_zero_loc, upper_zero_loc))

        # Split the lobe in 9 lobe points
        x = np.linspace(lower_zero_loc * 1.001, upper_zero_loc * 0.999, 9)

        # Store the relevant lobe
        f_interp = interp1d(range(nfft), mod_fft_s)
        mag_lobe = f_interp(x)
        f_interp = interp1d(lobe_index, np.unwrap(np.angle(fft_s[lobe_index])))
        arg_lobe = f_interp(x)
        x_lobe = ((x - 1) / nfft) - 0.0227  # f0/fs = 0.0227