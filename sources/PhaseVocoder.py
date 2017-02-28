#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

# Thi file contains the definition of the Phase vocoder class

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
    """The abstract class PhaseVocoder is used to get a phase-vocoder spectrum."""
    def __init__(self, analysis_hop, synthesis_hop, current_analysis_spectrum):
        if analysis_hop >= current_analysis_spectrum.get_nfft() or \
                synthesis_hop >= current_analysis_spectrum.get_nfft():
            raise InconsistentHopSizeError('Synthesis and/or Analysis Hop provided have inconsistent size')
        else:
            self._analysis_hop = analysis_hop
            self._synthesis_hop = synthesis_hop
            self._omega = np.array(range(current_analysis_spectrum.get_nfft())) *\
                ((2*np.pi)/current_analysis_spectrum.get_nfft())
            self._past_analysis_spectrum = Spectrum.void_spectrum(current_analysis_spectrum.get_nfft())
            self._past_synthesis_spectrum = Spectrum.void_spectrum(current_analysis_spectrum.get_nfft())
            self._current_analysis_spectrum = Spectrum.void_spectrum(current_analysis_spectrum.get_nfft())
            self.current_synthesis_spectrum = Spectrum.void_spectrum(current_analysis_spectrum.get_nfft())

    def _set_current_analysis_spectrum(self, new_spectrum):
        print "Set phase vocoder"
        # Assumed correct scheme
        self._past_analysis_spectrum = self._current_analysis_spectrum
        self._past_synthesis_spectrum = self.current_synthesis_spectrum
        self.current_synthesis_spectrum = Spectrum.void_spectrum(self._current_analysis_spectrum.get_nfft())
        self._current_analysis_spectrum = new_spectrum
        # Test scheme
        # self._past_analysis_spectrum = self._current_analysis_spectrum
        # self._past_synthesis_spectrum = new_spectrum
        # self._current_analysis_spectrum = self.current_synthesis_spectrum
        # self.current_synthesis_spectrum = Spectrum.void_spectrum(self._current_analysis_spectrum.get_nfft())
    current_analysis_spectrum = property(fset=_set_current_analysis_spectrum)

    @abstractmethod
    def get_pv_spectrum(self):
        pass

    def _kronecker_array(self, k):
        ek = np.zeros(self._current_analysis_spectrum.get_nfft())
        ek[k] = 1
        return ek

# Stationary Phase Vocoder


class StationaryPhaseVocoder(PhaseVocoder):

    def get_pv_spectrum(self):
        """Phase vocoder algorithm"""

        for k in xrange(self.current_synthesis_spectrum.get_nfft()):
            amplitude = self._current_analysis_spectrum.get_amplitude(k)  # \
                # if self._current_analysis_spectrum.get_amplitude(k) > 1e-12 else 1e-12
            phase = self._current_analysis_spectrum.get_phase(k)
            past_synth_phase = self._past_synthesis_spectrum.get_phase(k)
            past_analysis_phase = self._past_analysis_spectrum.get_phase(k)
            nfft = self._past_analysis_spectrum.get_nfft()

            # Get the phase difference
            delta_phi = phase - past_analysis_phase

            # Remove the expected phase difference
            # Note : (2 * np.pi * k) / nfft = omega(k)
            # delta_phi_prime = delta_phi - self._analysis_hop * (2 * np.pi * k) / nfft
            delta_phi_prime = delta_phi - self._analysis_hop * self._omega[k]

            # Map to - pi / pi range
            delta_phi_prime_mod = np.mod((delta_phi_prime + np.pi), (2 * np.pi)) - np.pi

            # Get the true frequency
            # true_freq = (2 * np.pi * k) / nfft + delta_phi_prime_mod / self._analysis_hop
            true_freq = self._omega[k] + delta_phi_prime_mod / self._analysis_hop

            # Get the final phase
            self.current_synthesis_spectrum.set_spectrum(np.array([amplitude]),
                                                         np.array([past_synth_phase + self._synthesis_hop * true_freq]),
                                                         start_bin=k, stop_bin=k+1)
        return self.current_synthesis_spectrum


class NonStationaryPhaseVocoderScalePhaseLocking(PhaseVocoder):

    def get_pv_spectrum(self):
        """Phase vocoder algorithm : Scale Phase-Locking"""

        for k in xrange(self.current_synthesis_spectrum.get_nfft()): #todo : peak finder + peak trajectory
            # todo : get_regions()

        # k1 is the current peak, k0 is the past peak

            amplitude = self._current_analysis_spectrum.get_amplitude(k)
            past_synth_phase = self._past_synthesis_spectrum.get_phase(k)
            past_analysis_phase_k0 = self._past_analysis_spectrum.get_phase(k)
            nfft = self._past_analysis_spectrum.get_nfft()

            # Get the phase difference
            delta_phi = phase - past_analysis_phase_k0

            # Remove the expected phase difference
            delta_phi_prime = delta_phi - self._analysis_hop * self._omega[k] #omega(k1)

            # Map to - pi / pi range
            delta_phi_prime_mod = (delta_phi_prime + np.pi) % (2 * np.pi) - np.pi

            # Get the true frequency
            true_freq = self._omega[k] + delta_phi_prime_mod / self._analysis_hop

            # Get phase for k1
            current_synthesis_phase_k1 = past_synth_phase + self._synthesis_hop * true_freq

            # Get phase for all k in the current frame region
            current_synthesis_phase = current_synthesis_phase_k1 + beta* (past_synth_phase - phase_k1)

            # Get the final phase in the region
            self.current_synthesis_spectrum += amplitude*np.exp(1j*current_synthesis_phase) *\
                self._kronecker_array(k)
        return self.current_synthesis_spectrum

        pass