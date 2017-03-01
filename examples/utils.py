#!/usr/bin/python
# -*- coding: utf-8 -*-

# Stationary and non-stationary sinusoidal model synthesis with phase vocoder and FFT-1
# Clément CAZORLA, Vincent CHRUN, Bastien FUNDARO, Clément MALIET

from scipy import signal
from scipy.interpolate import interp1d
import numpy as np
from sources.next2pow import next2pow
from sources.data_structure import *


class ParametersBuilder(object):
    def __init__(self, x, fs, analysis_hop, window_length, window_type, zero_padding_factor, threshold):
        self._x = x
        self._fs = fs
        self._analysis_hop = analysis_hop
        self._window_size = int(round(self._fs * window_length) if round(self._fs * window_length) % 2 != 0
                                else round(self._fs * window_length) + 1)
        self._window_type = window_type
        self._window = signal.get_window(self._window_type, self._window_size)
        self._window /= sum(self._window)
        self.number_frames = int(np.floor((len(self._x)-self._window_size)/analysis_hop))
        self._parameters_table = [None]*self.number_frames
        self._nfft = 2**(next2pow(self._window_size) + zero_padding_factor)
        self._t = 10**(threshold/20)
        self._gen_table()

    def _gen_table(self):
        frame = np.zeros(self._window_size)
        frame[:self._window_size] = self._x[:self._window_size]
        self._parameters_table[0] = self._get_lobes_param(frame, 0, phases=True)
        for i in xrange(1, self.number_frames):
            frame[:self._window_size] = self._x[i * self._analysis_hop:i * self._analysis_hop + self._window_size]
            self._parameters_table[i] = self._get_lobes_param(frame, i)

    def _get_lobes_param(self, frame, i, phases=False):
        frame *= self._window
        half_nfft = self._nfft / 2 + 1

        sw_gt = np.zeros(self._nfft)
        sw_gt[:(self._window_size - 1) / 2] = frame[((self._window_size + 1) / 2):]
        sw_gt[self._nfft - (self._window_size - 1) / 2 - 1:] = frame[:(self._window_size + 1) / 2]

        frame_spectrum = Spectrum.from_complex_spectrum(np.fft.fft(sw_gt, self._nfft))
        maximas = np.nonzero(np.logical_and(np.logical_and((frame_spectrum.amplitude[1:half_nfft-1] >
                                                            frame_spectrum.amplitude[:half_nfft-2]),
                             (frame_spectrum.amplitude[1:half_nfft-1] > frame_spectrum.amplitude[2:half_nfft])),
                             (frame_spectrum.amplitude[1:half_nfft - 1] > self._t)))[0]
        amplitudes = []
        frequencies = []
        phase = []
        for ploc in maximas:
            [iploc, ipmag, ipphase] = self._peak_interp(frame_spectrum, ploc)
            amplitudes.append(2*ipmag)
            frequencies.append(iploc/float(self._nfft))
            if phases:
                phase.append(ipphase)
            else:
                phase.append((i*self._analysis_hop*2*np.pi)*(iploc/float(self._nfft)))
        return Parameters(amplitudes, frequencies, phase)

    def _peak_interp(self, frame_spectrum, ploc):
        val = frame_spectrum.amplitude[ploc]
        lval = frame_spectrum.amplitude[ploc - 1]
        rval = frame_spectrum.amplitude[ploc + 1]
        iploc = ploc + .5 * ((lval - rval) / float(lval - 2 * val + rval))
        ipmag = val - .25 * ((lval - rval) * float(iploc - ploc))
        f = interp1d(np.arange(self._nfft), frame_spectrum.phase, "linear")
        ipphase = f(iploc) - 2*np.pi*(iploc/float(self._nfft))*self._analysis_hop
        return [iploc, ipmag, ipphase]

    def get_parameter(self, k):
        return self._parameters_table[k]
