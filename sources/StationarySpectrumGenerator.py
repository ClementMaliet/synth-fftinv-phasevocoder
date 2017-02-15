#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

# Import :

from SpectrumGenerator import *
from StationaryLobe import *


class StationarySpectrumGenerator(SpectrumGenerator):
    def __init__(self, window_type, window_size, parameters):
        SpectrumGenerator.__init__(self, window_type, window_size, parameters)
        self._lobe_generator = StationaryLobe(window_type, window_size)

    def _add_lobe(self, k):
        lobe = self._lobe_generator.get_lobe()
        amplitude = self._parameters.get_amplitude(k)
        phase = self._parameters.get_phase(k)
        nfft = self._spectrum.get_nfft()
        frequency = self._parameters.get_frequency(k) * nfft
        spectrum = Spectrum.void_spectrum(nfft)

        # size of positive freq. spectrum
        h_n = nfft / 2 + 1

        #    avoid frequencies out of range
        if (frequency > 0) or (frequency < h_n - 2):

            #   spectrum bins to fill
            b = np.int_(np.arange(9) - 4 + round(frequency) + 1)

            for m in xrange(9):
                #           peak lobe croses DC bin
                if b[m] < 0:
                    spectrum.set_spectrum(lobe.get_amplitude(m), phase, start_bin=-b[m], stop_bin=-b[m]+1)
                #           peak lobe croses Nyquist bin
                elif b[m] > (h_n - 1):
                    spectrum.set_spectrum(lobe.get_amplitude(m), phase, start_bin=2 * (h_n - 1)-b[m], stop_bin=2 *
                                          (h_n - 1) - b[m]+1)

                #           peak lobe in positive freq. range
                else:
                    spectrum.set_complex_spectrum(lobe.get_amplitude(m)*np.exp(1j*phase)+lobe.get_amplitude(m) *
                                                  np.exp(-1j*phase)*(b[m] == 0 or b[m] == h_n-1),
                                                  start_bin=b[m], stop_bin=b[m]+1)
            spectrum *= np.array(amplitude)
            self._spectrum += spectrum

        else:
            pass
