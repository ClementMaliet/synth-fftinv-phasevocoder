# -*- coding: utf8-*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

# Thi file contains the definition of the Phase vocoder class
import numpy as np
from abc import ABCMeta, abstractmethod
from data_structure import *

class PhaseVocoderError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass

class InconsistentHopSizeError(PhaseVocoderError):
    """Raised if the hop size is larger than the nfft size of given spectrum"""
    def __init__(self, message):
        self.message = message


class PhaseVocoder:
    __metaclass__ = ABCMeta
    """The abstract class PhaseVocoder is used to get a phase-vocoded spectrum."""
    def __init__(self, analysis_hop, synthesis_hop, omega, past_analysis_spectrum, past_synthesis_spectrum, current_analysis_spectrum, current_synthesis_spectrum):
        if len(analysis_hop) >= len(past_analysis_spectrum.get_nfft()) or len(synthesis_hop) >= len(past_synthesis_spectrum.get_nfft()):
            raise InconsistentHopSizeError('Synthesis and/or Analysis Hop provided have inconsistent size')
        else:
            self._analysis_hop = analysis_hop
            self._synthesis_hop = synthesis_hop
            self._omega = omega
            self._past_analysis_spectrum = past_analysis_spectrum
            self._past_synthesis_spectrum = past_synthesis_spectrum
            self._current_analysis_spectrum = current_analysis_spectrum
            self.current_synthesis_spectrum = current_synthesis_spectrum

    def _set_current_analysis_spectrum(self, new_spectrum):
        self._past_analysis_spectrum = self._current_analysis_spectrum
        self._past_synthesis_spectrum = self.current_synthesis_spectrum
        # current synth a zero
        self._current_analysis_spectrum = new_spectrum
    current_analysis_spectrum = property(_set_current_analysis_spectrum)

    @abstractmethod
    def get_pv_spectrum(self):
        pass


# Stationary Phase Vocoder

class StationaryPhaseVocoder(PhaseVocoder):
    def get_pv_spectrum(self):
        """Phase vocoder algorithm"""

        return current_synthesis_spectrum






"""
###############################
        # read input and get the timescale factor
        (sr, signalin) = wavfile.read(sys.argv[2])
        L = len(signalin)
        tscale = float(sys.argv[1])


        # signal blocks for processing and output
        phi = zeros(N)
        out = zeros(N, dtype=complex)
        sigout = zeros(L / tscale + N)


        # max input amp, window
        amp = max(signalin)
        win = hanning(N)
        p = 0
        pp = 0

        while p < L - (N + H):

        # take the spectra of two consecutive windows
        p1 = int(p)
        spec1 = fft(win * signalin[p1:p1 + N])
        spec2 = fft(win * signalin[p1 + H:p1 + N + H])


        # take their phase difference and integrate
        phi += (angle(spec2) - angle(spec1))


        # bring the phase back to between pi and -pi
        while phi < -pi: phi += 2 * pi
        while phi >= pi: phi -= 2 * pi
        out.real, out.imag = cos(phi), sin(phi)


        # inverse FFT and overlap-add
        sigout[pp:pp + N] += win * ifft(abs(spec2) * out)
        pp += H
        p += H * tscale """
