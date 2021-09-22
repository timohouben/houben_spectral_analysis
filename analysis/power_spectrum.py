# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------


def power_spectrum(input, output, time_step_size, method="scipyffthalf", o_i="oi"):
    """
    This script computes the power spectral density estimate of a time series.
    As default, output spectrum is devided by input spectrum.
    You can choose between three methods. 'scipyffthalf' and 'scipyperio'
    reveals almost exactly the same results. 'scipywelch' computes a smoothed
    periodogram.


    Parameters
    ----------
    input : 1D array, list
        Time series of and input process of e.g. a LTI system. If considering an
        aquifer as filter of the LTI, the input signal would be equal to the
        recharge time series of the aquifer.
    output : 1D array, list
        Time series of and input process of e.g. a LTI system. If considering an
        aquifer as filter of the LTI, the ouput signal would be euqual to the
        head time seires of the aquifer.
    time_step_size : integer
        The size of the time step between every data point in seconds.
    method : string, Default: 'scipyffthalf'
        Method which will be used to derive the spectrum.
        'scipyffthalf'
        # ======================================================================
        # method 1: Periodogram: Power Spectral Density: abs(X(w))^2
        #           http://staff.utia.cas.cz/barunik/files/QFII/04%20-%20Seminar/04-qf.html
        # ======================================================================
        'scipywelch'
        # ======================================================================
        # method 2: scipy.signal.welch
        #           https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.welch.html#r145
        # ======================================================================
        'scipyperio'
        # ======================================================================
        # method 3: Scipy.signal.periodogram
        #           https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.signal.periodogram.html
        # ======================================================================
    o_i : string
        'o_i' : output spectrum will be devided by input spectrum
        'i' : only input spectrum will be returned
        'o' : only output spectrum will be returned

    Yields
    ------
    frequency_xx : 1D array
        Corresponding frequencies of the Fourier Transform.

    power_spectrum_xx : 1D array
        Power spectrum of time series.

    Suggested Improvements
    ----------------------
    - Make input an optional argument! New structure: Every method is a function.


    """
    import numpy as np

    if np.shape(input) != np.shape(output) and o_i == "oi":
        raise ValueError("x and y must have same length.")
    if np.asarray(input).ndim != 1:
        raise ValueError("x and y must have dimension = 1.")
    len_input = len(input)
    len_output = len(output)
    # define the sampling frequency/time step
    # -------------------------------------------------------------------------
    sampling_frequency = 1.0 / time_step_size  # [Hz] second: 1, day: 1.1574074074074E-5
    # methodologies for power spectral density
    # -------------------------------------------------------------------------

    if method == "scipyffthalf_russian":
        import scipy.fftpack as fftpack
        # first value was popped because frequencies are very low (=0) and cause errors while fitting
        power_spectrum_input = fftpack.fft(input)
        power_spectrum_output = fftpack.fft(output)
        if len_input == len_output:
            power_spectrum_result = power_spectrum_output / power_spectrum_input
            power_spectrum_result = abs(power_spectrum_result[: int(round(len(power_spectrum_result) / 2))]) ** 2
            power_spectrum_result = power_spectrum_result[1:]
        frequency_input = (
            abs(fftpack.fftfreq(len_input, time_step_size))[
                : int(round(len_output / 2))
            ]
        )[1:]
        frequency_output = (
            abs(fftpack.fftfreq(len_output, time_step_size))[
                : int(round(len_output / 2))
            ]
        )[1:]

    elif method == "scipyffthalf":
        import scipy.fftpack as fftpack

        # first value was popped because frequencies are very low (=0) and cause errors while fitting
        spectrum = fftpack.fft(input)
        spectrum = abs(spectrum[: int(round(len(spectrum) / 2))]) ** 2
        power_spectrum_input = spectrum[1:]
        spectrum = fftpack.fft(output)
        spectrum = abs(spectrum[: int(round(len(spectrum) / 2))]) ** 2
        power_spectrum_output = spectrum[1:]
        if len_input == len_output:
            power_spectrum_result = power_spectrum_output / power_spectrum_input
        frequency_input = (
            abs(fftpack.fftfreq(len_input, time_step_size))[
                : int(round(len_input / 2))
            ]
        )[1:]
        frequency_output = (
            abs(fftpack.fftfreq(len_output, time_step_size))[
                : int(round(len_output / 2))
            ]
        )[1:]

    elif method == "scipywelch":
        from scipy import signal

        nperseg = int(round(len(input) / 10))
        frequency_input, power_spectrum_input = signal.welch(
            input, sampling_frequency, nperseg=nperseg, window="hamming"
        )
        frequency_output, power_spectrum_output = signal.welch(
            output, sampling_frequency, nperseg=nperseg, window="hamming"
        )
        if len_input == len_output:
            power_spectrum_result = power_spectrum_output / power_spectrum_input

    elif method == "scipyperio":
        from scipy import signal

        frequency_input, power_spectrum_input = signal.periodogram(
            input, fs=sampling_frequency
        )
        frequency_output, power_spectrum_output = signal.periodogram(
            output, fs=sampling_frequency
        )
        frequency_output = frequency_output[1:]
        frequency_input = frequency_input[1:]
        power_spectrum_input = power_spectrum_input[1:]
        power_spectrum_output = power_spectrum_output[1:]
        if len_input == len_output:
            power_spectrum_result = power_spectrum_output / power_spectrum_input
    else:
        print("Method not valid.")

    if o_i == "i":
        return np.asarray(frequency_input), np.asarray(power_spectrum_input)
        # return frequency_input, power_spectrum_input
    elif o_i == "o":
        return np.asarray(frequency_output), np.asarray(power_spectrum_output)
        # return frequency_input, power_spectrum_input
    elif o_i == "oi":
        return np.asarray(frequency_input), np.asarray(power_spectrum_result)
        # return frequency_input, power_spectrum_result
