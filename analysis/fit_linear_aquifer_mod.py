# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division
# ------------------------------------------------------------------------------

# Hey Estanis. The following script should do the job. It worked for me.
# Please try it and let me know if it works.
# Try to use python3 but with 2 it should work as well.

##############################################
# PLEASE SCROLL DOWN TO THE END OF THE FILE! #
##############################################


def discharge_ftf_linear_modified_lin(f, t_c, alpha):
    """
    This is the transfer function for the linear reservois model but slighly
    modified with an additional alpha to fit the slope of the spectra.

    Parameters
    ----------

    f : 1D-array, float
        The frequency, not angular frequency.
    t_c : scalar, float
        Characteristic time of the aquifer.
    alpha : float
        Added parameter to fit the slope of the spectra.
    """

    import numpy as np

    # calculate the angular frequency omega from f
    omega = [i * 2 * np.pi for i in f]

    spectrum = []
    for w in omega:
        spectrum.append((1 / ( 1 + w**2 * t_c**2))**alpha)
        #spectrum.append(1 / (1 + w**2 * t_c**2))
    spectrum = np.array(spectrum)
    return np.log10(spectrum)

def discharge_ftf_linear_modified_fit_lin(frequency, spectrum):
    """
    Function to fit the discharge_ftf_linear_modified to data.

    Parameters
    ----------

    frequency : 1D-array, float
        The frequency, not angular frequency.
    spectrum : 1D-array, float
        The PSD of the observation point.
    """
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(discharge_ftf_linear_modified_lin, frequency, spectrum)
    return popt, pcov


def discharge_ftf_linear_modified(f, t_c, alpha):
    """
    This is the transfer function for the linear reservois model but slighly
    modified with an additional alpha to fit the slope of the spectra.

    Parameters
    ----------

    f : 1D-array, float
        The frequency, not angular frequency.
    t_c : scalar, float
        Characteristic time of the aquifer.
    alpha : float
        Added parameter to fit the slope of the spectra.
    """

    import numpy as np

    # calculate the angular frequency omega from f
    omega = [i * 2 * np.pi for i in f]

    spectrum = []
    for w in omega:
        spectrum.append((1 / ( 1 + w**2 * t_c**2))**alpha)
        #spectrum.append(1 / (1 + w**2 * t_c**2))
    spectrum = np.array(spectrum)
    return spectrum

def discharge_ftf_linear_modified_fit(frequency, spectrum):
    """
    Function to fit the discharge_ftf_linear_modified to data.

    Parameters
    ----------

    frequency : 1D-array, float
        The frequency, not angular frequency.
    spectrum : 1D-array, float
        The PSD of the observation point.
    """
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(discharge_ftf_linear_modified, frequency, spectrum)
    return popt, pcov


# if __name__ == "__main__" just means, that it will be executed when you run the script.
if __name__ == "__main__":
    pass
    import numpy as np

    ##############################################
    # UNTIL HERE:                                #
    ##############################################

    # Please just change the following path:
    path_to_file_with_filename = "test_data_estanis.txt"

    data = np.loadtxt(path_to_file_with_filename)
    frequency = data[:,0]
    power = data[:,1]
    import matplotlib.pyplot as plt

    popt, pcov = discharge_ftf_linear_modified_fit(frequency,power)
    t_c = popt[0]
    print("t_c is : " + str(t_c))
    alpha = popt[1]
    print("alpha is: " + str(alpha))

    # fit on linear data
    popt_lin, pcov_lin = discharge_ftf_linear_modified_fit_lin(frequency,np.log10(power))
    t_c_lin = popt_lin[0]
    print("t_c_lin is : " + str(t_c_lin))
    alpha_lin = popt_lin[1]
    print("alpha_lin is: " + str(alpha_lin))

    plt.loglog(frequency, discharge_ftf_linear_modified(frequency, t_c, alpha), label="fit log")
    plt.loglog(frequency, np.power(10,discharge_ftf_linear_modified_lin(frequency, t_c_lin, alpha_lin)), label="fit lin")
    plt.loglog(frequency, power, label="data")
    plt.legend()
    plt.show()
