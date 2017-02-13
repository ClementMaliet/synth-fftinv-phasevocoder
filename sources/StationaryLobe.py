# Import the libraries

from scipy import signal
from data_structure import *
import numpy as np
from data_structure import *


def nextpow(i):
    n = 1
    while n < i: n *= 2
    return n
# Classes:

class StationaryLobe:
    def __init__(self, window_type, window_size):
        self._window_type = window_type
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._lobe = Spectrum([],[])
        self._gen_lobe()

    def _set_window_size(self, window_size):
        self._window_size = window_size
        self._window = signal.get_window(self._window_type, self._window_size)
        self._gen_lobe()

    def _set_window_type(self, window_type):
        self._window_type = window_type
        self._window = signal.get_window(self._window_type, self._window_size)
        self._gen_lobe()



    window_type = property(_set_window_type)
    window_size = property(_set_window_size)

    def _gen_lobe(self):
        w1 = np.fft.fft(self._window, self._window_size)
        w1 = np.concatenate([w1[-4:], w1[:5]])
        print w1
        self._lobe.set_complex_spectrum(w1)



    def get_lobe(self):
        return self._lobe
