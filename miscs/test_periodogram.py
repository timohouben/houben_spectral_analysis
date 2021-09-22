import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt
from scipy.signal import periodogram
import sys
sys.path.append("/Users/houben/phd/python/scripts/spectral_analysis/")
from analysis.calc_tc import calc_tc

signal1 = np.random.rand(100)


xvalue = np.linspace(1,100,1000)
sine = np.sin(xvalue)

plt.plot(sine)

signal = signal1

# via autocov
signal_acorr = np.correlate(signal,signal, "full")
plt.plot(signal_acorr)
spectrum1 = fft.fft(signal_acorr)
plt.plot(spectrum1)

# via magnitude squared
spectrum2 = fft.fft(signal)
spectrum2 = abs(spectrum2) ** 2
plt.plot(spectrum2)

# periodogram
f, Pxx = periodogram(signal)
plt.plot(Pxx)

plt.psd(signal)
plt.plot(spectrum2[1:])


### some calculations which has been used for paper 

S = 0.1
T = 0.01
L = 1000
tc = calc_tc(L, S, T, which="liang2013")

S = 0.01
T = 0.001
L = 5000
tc = calc_tc(L, S, T, which="liang2013")
