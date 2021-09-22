# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------
import numpy as np

def varq_analytical_sum_helper(m, n, plot=False):
    """
    Function to calculate the sum of varq_analytical for large numbers
    of n and m.
    """
    def b_n(n):
        return (2 * n + 1) * np.pi / 2

    def b_m(m):
        return (2 * m + 1) * np.pi / 2

    sum_evo = []

    outer_sum = 0
    for j in range(0, m):
        inner_sum = 0
        for k in range(0, n):
            inner_sum += ((-1) ** (j + k) * 8) / (b_n(k) ** 2 + b_m(j) ** 2)
        outer_sum += inner_sum
        sum_evo.append(outer_sum)
        if j % 500 == 0:
            print(j)

    if plot is True:
        import matplotlib.pyplot as plt
        plt.plot(range(len(sum_evo)), sum_evo)

    return outer_sum

def varq_analytical(sigma_n, T, S, L, m=None, n=None, tau_n=1):
    """
    Analytical solution for the baseflow variance.
    """
    # function variables
    alpha = T / S / L ** 2

    def b_n(n):
        return (2 * n + 1) * np.pi / 2

    def b_m(m):
        return (2 * m + 1) * np.pi / 2

    if m is not None and n is not None:
        # baseflow variance for provided m and n
        outer_sum = 0
        for j in range(0, m):
            inner_sum = 0
            for k in range(0, n):
                inner_sum += ((-1) ** (j + k) * 8) / (b_n(k) ** 2 + b_m(j) ** 2)
            outer_sum += inner_sum
        varq = sigma_n * tau_n * alpha * outer_sum
        return varq
    else:
        # value of varq_analytical_sum_helper for m, n = 20000
        #print("m and n == None! Set m and n to 200000")
        outer_sum = 1.2223146765429733
        varq = sigma_n * tau_n * alpha * outer_sum
        return varq



def varq_analytical_fit(var, sigma_n, T, S, L, m, n, tau_n_bounds=[1,10000000000], step=10000000):
    """
    Use this function to fit the analytical solution of the baseflow variance
    to a certain value (var) and optimize tau_n.

    Parameters
    ----------

    var : float
        Variance to get fitted to.
    sigma_n : float
        Variance of recharge.
        ...


    Return
    ------

    min_array : array of shape (1, 4)
        [0, 0] : optimal tau_n
        [0, 1] : variance of series to be fitted to
        [0, 2] : variance of analytical solution
        [0, 3] : difference between variance of series and variance of anal. sol.
    """

    min_array = []
    tau_ns = np.arange(tau_n_bounds[0], tau_n_bounds[1], step)
    for tau_n in tau_ns:
        #print("Testing for tau_n = " + str(tau_n))
        var_temp = varq_analytical(sigma_n, T, S, L, m, n, tau_n)
        min_array.append([tau_n, var, var_temp, abs(var-var_temp)])

    min_array = np.array(min_array)
    min_idx = np.where(min_array[:,3] == np.min(min_array[:,3]))
    #print(min_array)
    return min_array[min_idx,:]

def variance_evolution(series):
    """
    Return the evolution of the variance along the time series.

    Parameters
    ----------

    series : 1D array
        Series to calculate the variance from

    Returns
    -------

    variance : 1D array
        The variance evolution.
    varnorm : 1D array
        The normalized variance evolution. Max(varnorm) == 1
    """

    variance = []
    for timestep in np.arange(1, len(series) + 1):
        variance.append(np.var(series[0:timestep]))
    varnorm = [i / np.nanmax(variance) for i in variance]
    variance = np.array(variance)
    varnorm = np.array(varnorm)
    return variance, varnorm


if __name__ == "__main__":
    pass
    import matplotlib.pyplot as plt

    # test the function
    sigma_n = 0.01
    m = 1000
    n = 1000
    T = 0.01
    S = 0.1
    L = 1000
    time_step_size = 86400
    varq = []
    # for n_, m_ in zip(range(n), range(n)):
    #     # print("Currently calculating n/m: " + str(n_) + " " + str(m_))
    #     varq.append(varq_analytical(sigma_n, T, S, L, m_, n_, time_step_size))
    # plt.ylabel("variance")
    # plt.xlabel("iteration n, m")
    # plt.title("Convergence of Semi-Analytical Solution of Baseflow Variance")
    # plt.plot(varq)
    # plt.show()
    #
    s = varq_analytical(sigma_n, T, S, L, m, n, time_step_size)

    a = varq_analytical_fit(0.00010560, sigma_n, T, S, L, m, n, tau_n_bounds=[0,1], step=0.11)

    # coded in Wolfram Alpha
    # 85*86400*(0.01/0.1/1000^2)*sum (-1)^(n*m) * 8 / (((2*n+1)*pi/2)^2 + ((2*m+1)*pi/2)^2) from n=0 to 250, m=0 to 250
