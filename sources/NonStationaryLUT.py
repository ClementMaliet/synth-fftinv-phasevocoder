#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET


from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from data_structure import *
from LobeGenerator import LobeGenerator
from NonStationaryLobe import *


class NonStationaryLUT(LobeGenerator):
    def __init__(self, regular_grid, acr_domain, fcr_domain, number_acr, number_fcr, window_type, window_size, nfft,
                 fs=None, method=None):
        if fs is None:
            self._fs = 44100
        else:
            self._fs = fs
        if method is None:
            self._method = "linear"
        else:
            self._method = method
        LobeGenerator.__init__(self, window_type, window_size, nfft)
        self._abscisse = []
        self._ordonnee = []
        self._interpolated_lobe = NonStationaryLobe.void_lobe(9)
        self._regular_grid = regular_grid
        self._domain = [acr_domain, fcr_domain]
        self._number_acr = number_acr
        self._number_fcr = number_fcr
        self._number_points = number_acr*number_fcr
        self._lut = []
        self._gen_lobe()


    def _gen_uniform_lut(self):
        acr = np.linspace(self._domain[0][0], self._domain[0][1], self._number_acr)
        fcr = np.linspace(self._domain[1][0], self._domain[1][1], self._number_fcr)
        n = np.array(range(self._window_size))
        t = (1/self._fs)*n

        for i in xrange(self._number_acr):
            for j in xrange(self._number_fcr):
                lobe = self._gen_lobes_legacy(i, j, acr, fcr, t, n)
                self._abscisse.append(acr[i])
                self._ordonnee.append(fcr[j])
                self._lut.append(lobe)

    def _gen_non_uniform_lut(self):
        pass

    def _gen_lobes_legacy(self, i, j, acr, fcr, t, n):
        # Generate chirp
        mu = acr[i]
        psi = fcr[j]
        s = (1 + mu * t) * np.exp(1j * ((0.5 * psi * t ** 2) + (2 * np.pi * 0.0227 * n)))

        # Apply window
        s *= self._window

        # Correct phase
        sw = np.zeros(self._nfft)
        sw[:(self._window_size - 1) / 2] = s[(self._window_size + 1) / 2:]
        sw[self._nfft - (self._window_size - 1) / 2:] = s[:(self._window_size + 1) / 2 - 1]

        # Compute the fft
        fft_s = np.fft.fft(sw, self._nfft)
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
        f_interp = interp1d(range(self._nfft), mod_fft_s)
        amplitude = f_interp(x)
        f_interp = interp1d(lobe_index, np.unwrap(np.angle(fft_s[lobe_index])), fill_value = "extrapolate")
        phase = f_interp(x)
        abscisse = ((x - 1) / self._nfft) - 0.0227  # f0/fs = 0.0227
        lobe = NonStationaryLobe(amplitude, phase, abscisse)
        return lobe

    def _gen_lobe(self):
        if self._regular_grid:
            self._gen_uniform_lut()
        else:
            self._gen_non_uniform_lut()
        self._f = [interp2d(self._abscisse, self._ordonnee, [lobe.get_amplitude(k) for lobe in self._lut], self._method)
                   for k in range(9)]
        self._g = [interp2d(self._abscisse, self._ordonnee, [lobe.get_phase(k) for lobe in self._lut], self._method)
                   for k in range(9)]
        self._h = [interp2d(self._abscisse, self._ordonnee, [lobe.get_abscissa(k) for lobe in self._lut], self._method)
                   for k in range(9)]

    def get_lobe(self):
        return self._interpolated_lobe


    def interpolate_lobe(self, acr, fcr, method=None):
        if method != self._method:

            self._method = method
            self._f = [
                interp2d(self._abscisse, self._ordonnee, [lobe.get_amplitude(k) for lobe in self._lut], self._method)
                for k in range(9)]
            self._g = [interp2d(self._abscisse, self._ordonnee, [lobe.get_phase(k) for lobe in self._lut], self._method)
                       for k in range(9)]
            self._h = [
                interp2d(self._abscisse, self._ordonnee, [lobe.get_abscissa(k) for lobe in self._lut], self._method)
                for k in range(9)]

        lobe_mag = np.array([fk(acr, fcr) for fk in self._f])
        lobe_phase = np.array([gk(acr, fcr) for gk in self._g])
        lobe_frequency = np.array([hk(acr, fcr) for hk in self._h])

        self._interpolated_lobe = NonStationaryLobe(lobe_mag, lobe_phase, lobe_frequency)


