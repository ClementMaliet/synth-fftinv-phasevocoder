#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from SpectrumGenerator import *
from NonStationaryLobeGenerator import *  # todo : LobeGeneration in NonStationaryLobeGenerator


class NonStationarySpectrumGeneratorError(SpectrumGeneratorError):
    """Base class for exception regarding the non-stationary spectrum generator class"""
    pass


class NonStationarySpectrumGenerator(SpectrumGenerator):
    def __init__(self, window_type, window_size, parameters, nfft, analysis_hop, regular_grid,
                 acr_domain, fcr_domain, number_acr, number_fcr):
        self._regular_lut = regular_grid
        self._lobe_generator = NonStationaryLobeGenerator(regular_grid, acr_domain, fcr_domain, number_acr,
                                                          number_fcr, window_type, window_size, nfft)
        SpectrumGenerator.__init__(self, window_size, parameters, nfft, analysis_hop)

    def _add_lobe(self, k, lobe):
        amplitude = self._parameters.get_amplitude(k) / 2.
        acr = self._parameters.get_acrs(k)
        fcr = self._parameters.get_fcrs(k)
        self._lobe_generator.interpolate_lobe(acr, fcr)
        lobe = self._lobe_generator.get_lobe()
        phase = self._parameters.get_phase(k) + 2*np.pi*self._parameters.get_frequency(k)*self._analysis_hop + \
            fcr * self._analysis_hop**2 / 2.
        nfft = self._spectrum.get_nfft()
        frequency = self._parameters.get_frequency(k) * nfft
        spectrum = Spectrum.void_spectrum(nfft)
        spectrum_under_dc_bin = Spectrum.void_spectrum(nfft)
        spectrum_over_nyquist_bin = Spectrum.void_spectrum(nfft)

        # size of positive freq. spectrum
        h_n = nfft / 2 + 1

        #    avoid frequencies out of range
        if (frequency > 0) or (frequency < h_n - 2):
            #   spectrum bins to fill
            freqs = lobe.x * nfft + frequency
            f = interp1d(freqs, lobe.amplitude, "linear")
            g = interp1d(freqs, lobe.phase, "linear")
            b = np.arange(np.int_(np.ceil(np.min(freqs))), np.int_(np.floor(np.max(freqs))))
            lobe_amplitudes = f(b)
            lobe_phases = g(b)
            normal = np.nonzero(np.logical_and(b >= 0, b <= h_n - 1))[0]
            under_dc_bin = np.nonzero(b < 0)[0]
            over_nyquist_bin = np.nonzero(b > (h_n - 1))[0]

            # Normal case
            spectrum.set_complex_spectrum(lobe_amplitudes[normal]*np.exp(1j*phase) +
                                          lobe_amplitudes[normal] * np.exp(-1j*phase) *
                                          (np.logical_or(b[normal] == 0, b[normal] == h_n-1)),
                                          start_bin=np.min(b[normal]), stop_bin=np.max(b[normal])+1)
            if len(under_dc_bin) > 0:
                # Peak crosses DC
                spectrum_under_dc_bin.set_spectrum(np.fliplr([lobe_amplitudes[under_dc_bin]])[0],
                                                   -phase,
                                                   start_bin=np.min(-b[under_dc_bin]),
                                                   stop_bin=np.max(-b[under_dc_bin])+1)
                spectrum += spectrum_under_dc_bin
            if len(over_nyquist_bin) > 0:
                # Peak crosses Nyquist bin
                spectrum_over_nyquist_bin.set_spectrum(np.fliplr([lobe_amplitudes[over_nyquist_bin]])[0],
                                                       -phase,
                                                       start_bin=2 * (h_n - 1) - np.max(b[over_nyquist_bin]),
                                                       stop_bin=2 * (h_n - 1) - np.min(b[over_nyquist_bin]) + 1)
                spectrum += spectrum_over_nyquist_bin

            spectrum *= np.array(amplitude)
            # We add the lobe to the spectrum
            self._spectrum += spectrum
            self._spectrum += np.concatenate((np.array([0]), np.conj(spectrum.get_complex_spectrum()[:0:-1])))

        else:
            pass