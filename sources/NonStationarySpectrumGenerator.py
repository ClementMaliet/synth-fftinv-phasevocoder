from SpectrumGenerator import *

class NonStationarySpectrumGeneratorError(Exception):
    """Base class for exception regarding the spectrum class"""
    pass

class NonStationarySpectrumGenerator(SpectrumGenerator):

    def _init_(self,_regular_lut,_lut,_window_size, parameters, spectrum):
        self._regular_lut = _regular_lut
        self._lut = _lut
        SpectrumGenerator.__init__(self,_window_size, parameters, spectrum)

    @classmethod
    def _add_lobe(self,k):
        pass