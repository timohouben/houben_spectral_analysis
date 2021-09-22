"""
Metrics to evaluate the goodness of fit of two spectra.
"""
import numpy as np

def rmse(series1, series2, trans=None):
    """
    Function to calculate the r2 value of two 1D series.

    Parameters
    ----------

    series1, series2 : 1D array-like
    trans : 'log'
        Data will be log10 transformed before calculating the error.

    Returns
    -------

    rmse : float
        Root Mean Squared Error of the given series.
    """

    if trans == 'log':
        return np.sqrt(((np.log10(series1) - np.log10(series2)) ** 2).mean())
    elif trans == None:
        return np.sqrt(((series1 - series2) ** 2).mean())
    else:
        raise ValueError("Invalid string given for keyword 'trans'.")


def rsquare(series1, series2, trans=None):
    """
    Function to calculate the r2 value of two 1D series.

    Parameters
    ----------

    series1, series2 : 1D array-like
    trans : 'log'
        Data will be log10 transformed before calculating the error.

    Returns
    -------

    r_squared : float
        R-Squared value of the input series.
    """

    if trans == 'log':
        series1 = np.log10(series1)
        series2 = np.log10(series2)
    elif trans == None:
        pass
    else:
        raise ValueError("Invalid string given for keyword 'trans'.")


    # calculate the r-squared value of the fit
    correlation_matrix = np.corrcoef(series1, series2)
    correlation_xy = correlation_matrix[0, 1]
    r_squared = correlation_xy ** 2

    return r_squared


if __name__ == "__main__":
    # load some sample data
    import os
    CWD = os.getcwd()
    CWD_p = os.path.dirname(CWD)
    spectrum_1 = np.loadtxt(os.path.join(CWD_p, "test_data", "well_36_spectrum.txt"))
    spectrum_2 = np.loadtxt(os.path.join(CWD_p, "test_data", "well_57a_spectrum.txt"))
    spectrum_3 = np.loadtxt(os.path.join(CWD_p, "test_data", "well_61_spectrum.txt"))
    spectrum_4 = np.loadtxt(os.path.join(CWD_p, "test_data", "well_74_spectrum.txt"))

    # spetrum 1 vs spectrum 2
    print("## spetrum 1 vs spectrum 2")
    rsquared_ = rsquare(spectrum_1, spectrum_2)
    print("R-squared is: ", rsquared_)

    rmse_ = rmse(spectrum_1, spectrum_2)
    print("rmse is: ", rmse_)
    
    rmse_ = rmse(spectrum_2, spectrum_1)
    print("rmse is: ", rmse_)

    import matplotlib.pyplot as plt
    plt.loglog(spectrum_1)
    plt.loglog(spectrum_2)
    plt.show()

    # logtransform the specrta to get a linear plot
    spectrum_1_log = np.log10(spectrum_1)
    spectrum_2_log = np.log10(spectrum_2)
    
    rsquared_ = rsquare(spectrum_1_log, spectrum_2_log)
    print("R-squared is: ", rsquared_)

    rmse_ = rmse(spectrum_1_log, spectrum_2_log)
    print("rmse is: ", rmse_)
    
    rmse_ = rmse(spectrum_2_log, spectrum_1_log)
    print("rmse is: ", rmse_)

    # spetrum 3 vs spectrum 4
    print("## spetrum 3 vs spectrum 4")
    rsquared_ = rsquare(spectrum_3, spectrum_4)
    print("R-squared is: ", rsquared_)

    rmse_ = rmse(spectrum_3, spectrum_4)
    print("rmse is: ", rmse_)
    
    rmse_ = rmse(spectrum_4, spectrum_3)
    print("rmse is: ", rmse_)

    plt.loglog(spectrum_3)
    plt.loglog(spectrum_4)
    plt.show()

    # logtransform the specrta to get a linear plot
    spectrum_3_log = np.log10(spectrum_3)
    spectrum_4_log = np.log10(spectrum_4)
    
    rsquared_ = rsquare(spectrum_3, spectrum_4, trans="log")
    print("R-squared is: ", rsquared_)

    rmse_ = rmse(spectrum_3_log, spectrum_4_log)
    print("rmse is: ", rmse_)
    
    rmse_ = rmse(spectrum_4, spectrum_3, trans='log')
    print("rmse is: ", rmse_)

    plt.plot(spectrum_3_log)
    plt.plot(spectrum_4_log)
    plt.show()

