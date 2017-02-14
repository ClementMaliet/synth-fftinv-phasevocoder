from SpectrumGenerator import *

class NonStationarySpectrumGenerator(Exception)
    """Base class for exception regarding the spectrum class"""
    pass

class NonStationarySpectrumGenerator:

    def _init_(self,_regular_lut,_lut):
        self._regular_lut = _regular_lut
        self._lut = _lut

    @classmethod
    def _add_lobe(self,k):
        pass