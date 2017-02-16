#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from abc import ABCMeta, abstractmethod
from NonStationaryLUT import *
from data_structure import *
from NonStationarySpectrumGenerator import *
from PhaseVocoder import *
from SpectrumGenerator import *
from StationaryLobe import *
from StationarySpectrumGenerator import *


def next2pow(x):
    return 2**int(np.ceil(np.log(float(x))/np.log(2.0)))


class Synthesizer(object):
    synth_number = 0
    __metaclass__ = ABCMeta

    def __init__(self, window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                 current_parameters, fs):
        assert isinstance(current_parameters, Parameters.__class__)
        if Synthesizer.synth_number >= 1:
            warnings.warn("More the one synthesizer will be instantiated, consider reconfiguring or resetting rather \
                          than starting from scratch", UserWarning)
        self._window_size = window_size
        self._window_type = window_type
        self._nfft = 2**(next2pow(window_size) + zero_padding_factor)
        self._analysis_hop = analysis_hop
        self._synthesis_hop = synthesis_hop
        self._current_parameters = current_parameters
        self._past_parameters = Parameters.void_parameters(current_parameters.get_number_sinuses())
        self._fs = fs
        self._past_spectrum = Spectrum.void_spectrum(self._nfft)
        self._current_spectrum = Spectrum.void_spectrum(self._nfft)
        self._spectrum_generator = SpectrumGenerator(window_type, window_size, current_parameters, self._nfft)
        self._phase_vocoder = PhaseVocoder(analysis_hop, synthesis_hop, self._current_spectrum)
        Synthesizer.synth_number += 1

    def __del__(self):
        Synthesizer.synth_number -= 1

    def reset_synthetizer(self, window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                 current_parameters, fs):
        assert isinstance(current_parameters, Parameters.__class__)
        self._window_size = window_size
        self._window_type = window_type
        self._nfft = 2**(next2pow(window_size) + zero_padding_factor)
        self._analysis_hop = analysis_hop
        self._synthesis_hop = synthesis_hop
        self._current_parameters = current_parameters
        self._past_parameters = Parameters.void_parameters(current_parameters.get_number_sinuses())
        self._fs = fs
        self._past_spectrum = Spectrum.void_spectrum(self._nfft)
        self._current_spectrum = Spectrum.void_spectrum(self._nfft)
        self._spectrum_generator = SpectrumGenerator(window_type, window_size, current_parameters, self._nfft)
        self._phase_vocoder = PhaseVocoder(analysis_hop, synthesis_hop, self._current_spectrum)

    @abstractmethod
    def get_next_frame(self):
        pass

    def set_next_frame(self, next_parameters):
        self._past_parameters = self._current_parameters
        self._past_spectrum = self._current_spectrum
        self._current_spectrum = Spectrum.void_spectrum(self._nfft)
        self._current_parameters = next_parameters

    def _set_window_size(self, window_size):
        self._window_size = window_size
        self._spectrum_generator.window_size = window_size
    window_size = property(fset=_set_window_size)

    def _set_window_type(self, window_type):
        self._window_type = window_type
        self._spectrum_generator.window_type = window_type
    window_type = property(fset=_set_window_type)

    def _set_analysis_hop(self, analysis_hop):
        self._analysis_hop = analysis_hop
        self._phase_vocoder._analysis_hop = analysis_hop
    analysis_hop = property(fset=_set_analysis_hop)

    def _set_synthesis_hop(self, synthesis_hop):
        self._synthesis_hop = synthesis_hop
        self._phase_vocoder._synthesis_hop = synthesis_hop
    synthesis_hop = property(fset=_set_synthesis_hop)


class StationarySynthesizer(Synthesizer):
    def __init__(self, window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                 current_parameters, fs):
        Synthesizer.__init__(self, window_size, window_type, zero_padding_factor, analysis_hop, synthesis_hop,
                             current_parameters, fs)
        self._spectrum_generator = StationarySpectrumGenerator(self._window_type, self._window_size,
                                                               self._current_parameters, self._nfft)
        self._phase_vocoder = StationaryPhaseVocoder(self._analysis_hop, self._synthesis_hop, self._current_spectrum)

    def get_next_frame(self):
        pass
