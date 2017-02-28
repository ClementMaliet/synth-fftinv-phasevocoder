# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 20:17:39 2015

@author: Emily Abia, Quentin Biache, Lili_ZHENG
This script contain functions used to perform parameters estimation of modulated
sinusoidal model with generalized derivative methode and reassignement methode.
"""

#==============================================================================
# import section
#==============================================================================
from numpy import *
from matplotlib.pyplot import *
from scipy import signal
import time



#==============================================================================
# function section
#==============================================================================

###############################################################################
#   parameters estimator of modulated sinusoidal model (derivative method)
#   reference: GENERALIZATION OF THE DERIVATIVE ANALYSIS METHOD TO NON-STATIONARY SINUSOIDAL MODELING
#              Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
#   input:
#       s: input signal
#       Fs: sampling frequency
#       w: analysis window
#       N_fft: fft size
#       th: threshold for peaks detection, in negative dB
#   output:
#       a1:amplitude
#       omega1:frequency in radian
#       phi1:initial phase
#       mu1:modulation of amplitude in s^-1
#       psi1:modulation of frequency in Hz/s
###############################################################################
def estimator1(s,Fs,w,N_fft,th):
#   calculate the first and second derivative of signal
    d1 = gradient(s)*float(Fs)
    d2 = gradient(d1)*float(Fs)

#   size of the signal and size of the fft
    N = size(s,0)

#   Time axis
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)

#   calculate the spectrum of the signal and it's first derivative
    X = fft.fft((w*s),N_fft)
    X1 = fft.fft((w*d1),N_fft)
    X2 = fft.fft((w*d2),N_fft)

    #   find the pic
    mX = 20*log10(abs(X[1:N_fft/2]))
    m = array(where((mX[1:-1]>mX[2:])& (mX[1:-1]>th) & (mX[1:-1]>mX[:-2])))+1
    m = m[0,:]

#   estimation of omega (frequency in radian)
    r1 = peak_value(X1,m)/peak_value(X,m)
    omega1 = imag(r1)

#   initiation of output array
    L = size(omega1)
    mu1 = zeros(L)
    psi1 = zeros(L)
    phi1 = zeros(L)
    a1 = zeros(L)

#perform estimation
    for i in xrange(L):
        X =sum(w*s*exp(-1j*omega1[i]*t))
        X1 =sum(w*d1*exp(-1j*omega1[i]*t))
        X2 =sum(w*d2*exp(-1j*omega1[i]*t))

        r1 = X1 / X
        r2 = X2 / X

        mu1[i] = real(r1)
        psi1[i] = imag(r2) - 2.0*mu1[i]*omega1[i]
        p = Gamma_est(N, Fs, 0, mu1[i], psi1[i], w)
        phi1[i] = angle(X/p)
        a1[i] = 2.0*abs(X/p)
    return [a1,omega1,phi1,mu1,psi1]






###############################################################################
#   parameters estimator of mudulated sinusoidal model (reassignement methode)
#   reference: Techniques for the Automated Analysis of Musical Audio
#              Stephen Webley Hainsworth, A thesis submitted to the University
#              of Cambridge
#   input:
#       s: input signal
#       Fs: sampling frequency
#       th: threshold for peaks detection, in negative dB
#   output:
#       a1:amplitude
#       omega1:frequency in radian
#       phi1:initial phase
#       mu1:modulation of amplitude in s^-1
#       psi1:modulation of frequency in Hz/s
###############################################################################
def estimator2(s,Fs,th):

#   time axis
    N = size(s)
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)

#   used window and its first and second derivative
    w = signal.blackmanharris(N)
    wd = win_d1(N,Fs)
    wdd = win_d2(N,Fs)
    tw = (t)*w
    twd = (t)*wd

#   compute the spectre
    Xw = fft.fft(w*s)
    Xwd = fft.fft(wd*s)
    Xtw = fft.fft(tw*s)
    Xwdd = fft.fft(wdd*s)
    Xtwd = fft.fft(twd*s)

#   find the pic
    mX = 20*log10(abs(Xw[1:hN/2]))
    m = array(where((mX[1:-1]>mX[2:])& (mX[1:-1]>th) & (mX[1:-1]>mX[:-2])))+1
    m = m[0,:]

#   initiation of output array
    L = size(m)
    delta_omega = zeros(L)
    base_omega = zeros(L)
    omega1 = zeros(L)
    mu1 = zeros(L)
    psi1 = zeros(L)
    phi1 = zeros(L)
    a1 = zeros(L)

#   perform estimation
    delta_omega = -imag (Xwd[m] / Xw[m])
    base_omega = (m)*2*pi*Fs/N;
    omega1[:] = base_omega + delta_omega
    mu1[:] = -real (Xwd[m] / Xw[m])
    psi1[:] = (imag(Xwdd[m]/Xw[m]) - imag((Xwd[m]/Xw[m])**2)) / (real((Xtw[m]*Xwd[m])/(Xw[m])**2) - real(Xtwd[m]/Xw[m]))
    p = Gamma_est(N, Fs, delta_omega, mu1, psi1, w)
    phi1[:] = angle (Xw[m] / p)
    a1[:] = 2*abs (Xw[m] / p)
    a1,omega1,phi1,mu1,psi1 = array(a1),array(omega1),array(phi1),array(mu1),array(psi1)
    return [a1,omega1,phi1,mu1,psi1]











###############################################################################
#   compute the first and second derivative of blackman harris window
#   N: window size
#   Fs: sampling frequency
###############################################################################
def win_d1(N,Fs):
#   time axis
    T = N/float(Fs)
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)

#   coef of blackman_harris window
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;

    wd = -a1*(2*pi/T)*sin(2*pi*t/T)-a2*(4*pi/T)*sin(4*pi*t/T)-a3*(6*pi/T)*sin(6*pi*t/T)

    return wd

def win_d2(N,Fs):
#   time axis
    T = N/float(Fs)
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)

#   coef of blackman_harris window
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;

    wdd = -a1*(2*pi/T)**2*cos(2*pi*t/T)-a2*(4*pi/T)**2*cos(4*pi*t/T)-a3*(6*pi/T)**2*cos(6*pi*t/T)

    return wdd











###############################################################################
#   interpolation of spectral peak
#   mX: magnitude spectrum, ploc: locations of peaks
#   peak_value: interpolated values
#   note that ploc values are assumed to be between 2 and length(mX)-1
###############################################################################
def peak_value(mX,ploc):

#   magnitude of peak bin
    val = mX[ploc]

#   magnitude of bin at lef
    lval = mX[ploc-1]

#    magnitude of bin at right
    rval = mX[ploc+1]

#   center of parapola
    iploc = ploc+0.5*(lval-rval)/(lval-2*val+rval)

#   magitude of peaks
    peak_value = val-0.25*(lval-rval)*(iploc-ploc)

    return peak_value











###############################################################################
#   Gamma function used to perform estimation
#   reference: GENERALIZATION OF THE DERIVATIVE ANALYSIS METHOD TO NON-STATIONARY SINUSOIDAL MODELING
#              Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
###############################################################################
def Gamma_est(N,Fs,Delta,mu,psi,w):
#   fit the size of each term
    L = size(Delta)
    Mw = tile(w,[L,1])
    mu = transpose(tile(mu,[N,1]))
    Delta = transpose(tile(Delta,[N,1]))
    psi = transpose(tile(psi,[N,1]))

#   time axis
    hN = (N-1)/2
    t = arange(-hN,hN+1)/float(Fs)
    t = tile(t,[L,1])
    result = transpose (sum (Mw*exp(mu*t+1j*(Delta*t+psi*(t**2)/2.0)),1))

#   in case of problems with large mu value...
    idx_nan = array(where(isnan(result)))
    idx_inf = array(where(isinf(result)))
    idx_nan,idx_inf=idx_nan[0,:],idx_inf[0,:]
    idx = concatenate((idx_nan,idx_inf))

#   reset to default: no correction
    result[idx] = 1
    return result






##############################################################################
#   Compute a spectrum from a series of sine values
#   iploc, ipmag, ipphase: sine locations, magnitudes and phases
#   N: size of complex spectrum
#   Y: generated complex spectrum of sines
###############################################################################
def genspecsines(ploc, pmag, pphase, N):

# initialize output spectrum
    Y = zeros(N,dtype=complex)

# size of positive freq. spectrum
    hN = N/2+1

# generate all sines lobes
    for i in xrange(size(ploc)):

#    location of peak
        loc = ploc[i]

#    avoid frequencies out of range
        if ((loc<=0) or (loc>=hN-2)):continue

        binremainder = round(loc) - loc

#   main lobe bin to read
        lb = arange(9)-4+binremainder

#   lobe magnitude of the complex exponential
        lmag = genbh92lobe(lb,N,N)*(10**(pmag[i]/20.0))

#   spectrum bins to fill
        b = int_(arange(9)-4+round(loc)+1)

        for m in xrange(9):
#           peak lobe croses DC bin
            if (b[m]<0):Y[-b[m]] += lmag[m]*exp(-1j*pphase[i])

#           peak lobe croses Nyquist bin
            elif (b[m]>hN-1):Y[2*(hN-1)-b[m]] +=  lmag[m]*exp(-1j*pphase[i])

#           peak lobe in positive freq. range
            else:
                Y[b[m]]+=(lmag[m]*exp(1j*pphase[i])+lmag[m]*exp(-1j*pphase[i])*(b[m]==0 or b[m]==hN-1))

#   fill the rest of the spectrum
    Y[hN:] = conjugate(Y[hN-2:0:-1])
    return Y










###############################################################################
# x: bin positions to compute (real values)
# y: transform values
# M: Window size
# N: fft size
###############################################################################
def genbh92lobe(x,M,N):

    x = array(x)

#   frequency sampling
    f = 2.0*x*pi/N
    df = 2.0*pi/N
    r = float(N)/M

#   initialize window
    y = zeros(size(x))

#   window constants
    consts = array([0.35875, 0.48829, 0.14128, 0.01168])

#    sum Dirichlet kernels
#    for m in xrange(4):
#        y = y + 0.5*consts[m]*(D(f-df*m,N)+D(f+df*m,N))
    y = 0.5*consts[0]*(D(f-df*0*r,M)+D(f+df*0*r,M)) \
      + 0.5*consts[1]*(D(f-df*1*r,M)+D(f+df*1*r,M)) \
      + 0.5*consts[2]*(D(f-df*2*r,M)+D(f+df*2*r,M)) \
      + 0.5*consts[3]*(D(f-df*3*r,M)+D(f+df*3*r,M))

#   normalize
    y = (y)/((M)*consts[0])
    return abs(y)

def D(x,M):
#   Calculate rectangular window transform (Dirichlet kernel)
    y = sin(M*x/2.0)/sin(x/2.0)
#    y = y*exp(-1j*x*(N-1)/2.0)

#   avoid NaN if x==0
    y[where(y!=y)]=M
    return y











#==============================================================================
#  Time consume
#==============================================================================
def tic():
    #Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"