# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------


def calc_tc(L, S, T, which="linear"):
    """
    Calculates tc characteristic time scale (Gelhar 1974) in days!

    Parameters
    ----------

    L : float
        Aquifer length [m]
    S : float
        Storativity [-]
    T : float
        Transmissivity [m^2/s]
    which : string
        "linear" : for linear aquifer (default)
        "dupuit" : for dupuit aquifer based on Liang and Zhang 2013
        "lian2015" : for dupuit aquifer based on Liang and Zhang 2015
        "gelhar" : fir linear aquifer based on Gelhar 1993 (TBC)
        "liang2013" : for dupuit aquifer based on Liang and Zhang 2013

    Yields
    ------

    tc : float
        characteristic time scale [day]
    """
    if which == "linear_alt":
        return L ** 2 * S / T / 86400
    if which == "linear" or which == "gelhar":
        return L ** 2 * S / 3 / T / 86400
    if which == "liang2015":
        return L ** 2 * S / T / 86400
    if which == "dupuit" or which == "liang2013":
        from numpy import pi
        return 4 * L ** 2 * S / pi ** 2 / T / 86400
    else:
        print("Parameter 'which' can only be 'linear', 'gelhar', 'liang2013', 'dupuit' or 'liang2015'.")


if __name__ == "__main__":

    T = 3e-4
    T = 3e-2
    S = 0.0001*30
    L = 1000
    a = calc_tc(L, S, T, which="dupuit")
    print(a)