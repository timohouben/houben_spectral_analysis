# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

"""
We decided to change the manually added samples (L,x) for Birkach since the 
estimated characteristic time is too large and not realistic. Since the manually
added samples were added to the automatic samples and the SpA was performed
afterwards, we first need to perfom the SpA seperately with the new samples.
This will be added to a new results.csv file which can be taken to plot the
sensitivity plots for the paper.

for output see 
/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/20210708_edit_l_x_birkach
"""
# ------------------------------------------------------------------------------
# import modules
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import string


"""
OLD DESCRIPTION FROM THE ORIGINAL RUN:


Script to perform a sensitivity analysis for a fit of the shh_anlytical
solution to 4 groundwater head time series in the Main catchment. The target parameters are
T, S, tc, Rsquared and RMSE of the fit. The variables parameters are L and x.
Constraints: L must not be bigger than x!

L and x are drawn randomly from a gaussian?? distribution.

This script runs properly on EVE and can be also executet on a local machine with mpi4py installed.
Run the following command for 2 cores:

mpirun -n 2 python3 20200728_application_spectral_anal_sensitivity_new.py /Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivi 2 

"""
print(
    "###########################################################################################"
)
print("Time series will be detrended.")
print("The full spectrum will be evaluated.")
print(
    "###########################################################################################"
)
# import modules
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpi4py import MPI
import os
sys.path.append("/Users/houben/PhD/python/scripts/spectral_analysis")
#plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
# own modules
from analysis.power_spectrum import power_spectrum
from plot.plot_power_spectra import plot_spectrum
from analysis.shh_analytical import shh_analytical_fit, shh_analytical
from analysis.calc_tc import calc_tc
from analysis.processing import detrend
from analysis.goodness_of_fit import rmse, rsquare


#############################################
# Set parameters for script
# comment your analysis
comment = "editing_birkach_"
# set root path
# specify the path to the root containing a folder called "input" with input files in subfolders "gwm" and "recharge"

# output_path
output_path = os.path.abspath("/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/20210708_edit_l_x_birkach")

# Load head
birkach_h = os.path.abspath("/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/input/gwm/birkach_quadratic.txt")
# load recharge
birkach_r = os.path.abspath("/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/input/recharge/4416000.0_5522000.0_birkach_recharge_2019.txt")

names_h = ["Birkach"]
names_r = ["Birkach"]
paths_h = [birkach_h]
paths_r = [birkach_r]
# Further parameters for the fit.
m = None
n = None
norm = False
convergence = 0.01
time_step_size = 86400
# cut higher frequencies than cut_freq_higher
cut_freq_higher = 1e-4
# cut lower frequencies than cut_freq_lower
cut_freq_lower = 1e-9
# True: use only the overlap of time sereies of recharge and head
combine_df = True
detrend_ts = False
#############################################

#############################################
# Scripts begins
#############################################


# workaround to get script running on loca machine
rank = 0
slots = 1

# new tuple
own_list = (7509, 2997)
Ls, xs = [own_list[0]], [own_list[1]]
# ------------------------------------------------------------------------------------
# set path to results incl file name of results
path_to_results_df = os.path.join(
    output_path, comment + "results.csv"
)

# initialize a dataframe to save the results to
pd.set_option("precision", 10)
columns = [
    "name_h",
    "name_r",
    "L",
    "x",
    "S_out",
    "T_out",
    "tc_out",
    "cov",
    "time_step_size",
    "cut_freq_higher",
    "cut_freq_lower",
    "rsquared",
    "rsquared_log",
    "rmse",
    "rmse_log",
    "combine_df",
    "detrend_ts",
]
results = pd.DataFrame(columns=columns)

# loop over all paths and names with different aquifer lengths and positions
for name_h, path_h, name_r, path_r in zip(names_h, paths_h, names_r, paths_r):
    # Load the data
    head_df = pd.read_csv(path_h, sep=" ", header=None, names=["date", "head"])
    head_df["date"] = pd.to_datetime(head_df["date"])

    # load recharge data
    recharge_df = pd.read_csv(path_r, sep=" ", header=None, names=["date", "recharge"])
    recharge_df["date"] = pd.to_datetime(recharge_df["date"])

    if combine_df is True:
        # combine dataframes and remove rows with nans
        combined_df = pd.merge_ordered(recharge_df, head_df, how="inner")
        date_min = combined_df["date"].min()
        date_max = combined_df["date"].max()
        period = combined_df["date"].max() - combined_df["date"].min()
        print(
            "Start/end/length of series where head measurements and recharge overlap:\n"
            + str(date_min)
            + " .../... "
            + str(date_max)
            + " .../... "
            + str(period)
        )
        recharge_time_series = combined_df["recharge"].tolist()
        head_time_series = combined_df["head"].tolist()
    else:
        recharge_time_series = recharge_df["recharge"].tolist()
        head_time_series = head_df["head"].tolist()
        # modify the time series so that both have same length
        # assume: recharge is longest
        recharge_time_series = recharge_time_series[-len(head_time_series) :]

    if detrend_ts is True:
        head_time_series = detrend(head_time_series)
        recharge_time_series = detrend(recharge_time_series)
    else:
        pass

    # convert mm/d to recharge along the aquifer in m2/s
    recharge_time_series = [i / 86400 / 1000 for i in recharge_time_series]

    # calculate the power spectrum: Shh, output to FIT with analy solution only!
    frequency_output, Shh = power_spectrum(
        input=recharge_time_series,
        output=head_time_series,
        time_step_size=time_step_size,
        method="scipyffthalf",
        o_i="o",
    )
    frequency_input, Sww = power_spectrum(
        input=recharge_time_series,
        output=head_time_series,
        time_step_size=time_step_size,
        method="scipyffthalf",
        o_i="i",
    )

    if rank == 0:
        print("Currently plotting spectrum of " + str(name_h) + " ...")
        # plot only spectrum of Shh
        plot_spectrum(
            [Shh],
            frequency_output,
            heading="Shh - Head Power Spectrum " + name_h,
            labels=["Shh obs"],
            path=output_path,
            linestyle=["-"],
            marker=[""],
            lims=[(2e-9, 7e-6), (1e-5, 1e7)],
            name="Shh_" + name_h,
            ext=".eps",
        )

        # plot only spectrum of Sww
        plot_spectrum(
            [Sww],
            frequency_input,
            heading="Sww - Recharge Power Spectrum  " + name_r,
            labels=["Sww mHM"],
            path=output_path,
            linestyle=["-"],
            marker=[""],
            lims=[(2e-9, 7e-6), (1e-20, 1e-9)],
            name="Sww_" + name_r,
            ext=".eps",
        )

    # cut higher frequencies than cut_freq_higher
    cut_array_higher = np.less(frequency_input, cut_freq_higher)
    Sww = Sww[cut_array_higher]
    Shh = Shh[cut_array_higher]
    frequency_input = frequency_input[cut_array_higher]
    frequency_output = frequency_output[cut_array_higher]
    # cut lower frequencies than cut_freq_lower
    cut_array_lower = np.invert(np.less(frequency_input, cut_freq_lower))
    Sww = Sww[cut_array_lower]
    Shh = Shh[cut_array_lower]
    frequency_input = frequency_input[cut_array_lower]
    frequency_output = frequency_output[cut_array_lower]

    # fit the power spectrum with the analytical solution for every x and L
    for i, (L, x) in enumerate(zip(Ls, xs)):
        if i % slots == rank:
            print(
                "Currently fitting spectrum of "
                + str(name_h)
                + " for L = {:.0f} and x = {:.0f}".format(L, x)
                + " on rank "
                + str(rank)
            )
            try:
                popt, pcov = shh_analytical_fit(
                    Sww=Sww,
                    Shh=Shh,
                    f=frequency_input,
                    x=x,
                    m=m,
                    n=n,
                    L=L,
                    norm=False,
                    convergence=convergence,
                )
            except RuntimeError:
                print("Optimal parameters not found...")
                popt, pcov = [np.nan, np.nan], [[np.nan, np.nan], [np.nan, np.nan]]
                print("popt and pcov have been set to np.nan")
            except ValueError:
                print(
                    "either ydata or xdata contain NaNs, or if incompatible options are used"
                )
                popt, pcov = [np.nan, np.nan], [[np.nan, np.nan], [np.nan, np.nan]]
            except OptimizeWarning:
                print("Covariance of the parameters could not be estimated.")
                # popt, pcov = [np.nan, np.nan], [[np.nan, np.nan],[np.nan, np.nan]]

            Sy = abs(popt[0])
            T = abs(popt[1])

            # calculate the fitted power spectra
            Shh_fitted = shh_analytical(
                (frequency_input, Sww), Sy, T, x, L, m=n, n=m, norm=norm
            )

            # get the characteristic time
            tc = calc_tc(L, Sy, T, which="dupuit")
            # define a (discharge constant)
            # a = np.pi ** 2 * T / (4 * L ** 2)
            # define tc (characteristic time scale)
            # tc = Sy / a

            # calculate goodness of fit metrics
            rsquared = rsquare(Shh, Shh_fitted)
            rsquared_log = rsquare(Shh, Shh_fitted, trans="log")
            rmse_ = rmse(Shh, Shh_fitted)
            rmse_log = rmse(Shh, Shh_fitted, trans="log")

            data = np.vstack((Shh, Shh_fitted))

            labels = ["Shh numerical", "Shh fitted"]
            linestyle = ["-", "-"]
            marker = ["", "d"]
            figtxt = (
                "Derived Parameter:    S = %1.3e, T = %1.3e [m2/s], tc = %1.3e [d]\nInput Parameter:        L = %0.0f, x = %0.0f\nRsquared = %1.5f, Rsquared log = %1.5f, RMSE = %1.3e, RMSE log = %1.3f"
                % (Sy, T, tc, L, x, rsquared, rsquared_log, rmse_, rmse_log)
            )

            # plot Shh and the fitted spectrum
            plot_spectrum(
                data,
                frequency_input,
                heading="Shh - Head Power Spectrum " + name_h,
                labels=["Shh obs", "Shh fit"],
                path=output_path,
                linestyle=["-", " "],
                marker=["", "."],
                figtxt=figtxt,
                lims=[(2e-9, 7e-6), (1e-7, 1e8)],
                name="Shh_fitted_"
                + name_h
                + "_R_"
                + name_r
                + "_L_"
                + str(L)
                + "_x_"
                + str(x),
                ext=".eps",
            )
            
            # save the spectra
            np.savetxt(os.path.join(output_path, "Shh_" + name_h + ".txt"), Shh)
            np.savetxt(os.path.join(output_path, "Sww_" + name_r + ".txt"), Sww)
            np.savetxt(os.path.join(output_path, "Frequency_output_" + name_h + ".txt"), frequency_output)
            np.savetxt(os.path.join(output_path, "Frequency_input_" + name_r + ".txt"), frequency_input)
            np.savetxt(os.path.join(output_path, "Shh_fitted_"
                + name_h
                + "_R_"
                + name_r
                + "_L_"
                + str(L)
                + "_x_"
                + str(x)
                + ".txt"),
                Shh_fitted)

            results_temp = {
                "name_h": name_h,
                "name_r": name_r,
                "L": L,
                "x": x,
                "S_out": Sy,
                "T_out": T,
                "tc_out": tc,
                "cov": pcov,
                "time_step_size": time_step_size,
                "cut_freq_higher": cut_freq_higher,
                "cut_freq_lower": cut_freq_lower,
                "rsquared": rsquared,
                "rsquared_log": rsquared_log,
                "rmse": rmse_,
                "rmse_log": rmse_log,
                "combine_df": combine_df,
                "detrend_ts": detrend_ts,
            }
            results = results.append(other=results_temp, ignore_index=True, sort=False)
            # save temporary dataframe as file in between
            results.to_csv(path_to_results_df)
        else:
            continue
