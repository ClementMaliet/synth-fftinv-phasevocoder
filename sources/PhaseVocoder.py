#!/usr/bin/python
# -*- coding: utf-8 -*-

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


class PhaseVocoder(object):
    __metaclass__ = ABCMeta
    """The abstract class PhaseVocoder is used to get a phase-vocoded spectrum."""
    def __init__(self, analysis_hop, synthesis_hop, past_analysis_spectrum, past_synthesis_spectrum, current_analysis_spectrum, current_synthesis_spectrum):
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
        self._past_synthesis_spectrum = self._current_synthesis_spectrum
        # current synth a zero
        self._current_analysis_spectrum = new_spectrum
    current_analysis_spectrum = property(_set_current_analysis_spectrum)

    @abstractmethod
    def get_pv_spectrum(self):
        pass


# Stationary Phase Vocoder

class StationaryPhaseVocoder(PhaseVocoder):



    def get_pv_spectrum(self,k):
        """Phase vocoder algorithm"""
        amplitude = _current_analysis_spectrum.get_amplitude(k)
        phase = _current_analysis_spectrum.get_phase(k)
        previous_phase = _past_analysis_spectrum.get_phase(k)
        win_size = len(_past_analysis_spectrum.get_nfft(k))

        delta_phi = phase - previous_phase
        delta_phi_prime = delta_phi - self._analysis_hop * 2 * np.pi * np.array[0:(win_size - 1)]/win_size
        delta_phi_prime_mod = (delta_phi_prime + np.pi) % (2 * np.pi) - np.pi
        true_freq = 2 * np.pi * np.array[0:(win_size - 1)] / win_size + delta_phi_prime_mod / self._analysis_hop


        current_synthesis_spectrum *= phase + self._synthesis_hop * true_freq

        return current_synthesis_spectrum




