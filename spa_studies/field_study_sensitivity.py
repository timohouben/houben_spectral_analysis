"""
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

# set matplotlib style
CWD = os.getcwd()
CWD_p = os.path.dirname(CWD)
#plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
# add search path for own modules
sys.path.append("/Users/houben/PhD/python/scripts/spectral_analysis")
# add search path for own modules on eve
sys.path.append("/home/houben/python/scripts/spectral_analysis")
# own modules
from analysis.power_spectrum import power_spectrum
from plot.plot_power_spectra import plot_spectrum
from analysis.shh_analytical import shh_analytical_fit, shh_analytical
from analysis.calc_tc import calc_tc
from analysis.processing import detrend
from analysis.goodness_of_fit import rmse, rsquare

# some parameters for the mpi run
# get the number of slots from a system argument
slots = int(sys.argv[2])
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
#############################################
# Set parameters for script
# comment your analysis
comment = "1_"
# set root path
# specify the path to the root containing a folder called "input" with input files in subfolders "gwm" and "recharge"
try:
    root_path = sys.argv[1]
except IndexError:
    print("You forgot to give the path to multiple projects as argument...")
    root_path = input("Insert path to multiple projects: ")
# output_path
output_path = os.path.join(root_path, "output")
# make output directory
if rank == 0:
    if not os.path.exists(output_path):
        os.mkdir(output_path)
# Load head
birkach_h = os.path.join(root_path, "input/gwm/birkach_quadratic.txt")
stegaurach_h = os.path.join(root_path, "input/gwm/stegaurach_quadratic.txt")
strullendorf_west_h = os.path.join(root_path, "input/gwm/strull_west_quadratic.txt")
strullendorf_nord_h = os.path.join(root_path, "input/gwm/strull_nord_quadratic.txt")
# load recharge
birkach_r = os.path.join(
    root_path, "input/recharge/4416000.0_5522000.0_birkach_recharge_2019.txt"
)
stegaurach_r = os.path.join(
    root_path, "input/recharge/4416000.0_5526000.0_stegaurach_recharge_2019.txt"
)
strullendorf_west_r = os.path.join(
    root_path, "input/recharge/4424000.0_5522000.0_strullendorf_west_recharge_2019.txt"
)
strullendorf_nord_r = os.path.join(
    root_path, "input/recharge/4424000.0_5526000.0_strullendorf_nord_recharge_2019.txt"
)
names_h = ["Birkach", "Stegaurach", "Strullendorf_West", "Strullendorf_Nord"]
names_r = ["Birkach", "Stegaurach", "Strullendorf_West", "Strullendorf_Nord"]
paths_h = [birkach_h, stegaurach_h, strullendorf_west_h, strullendorf_nord_h]
paths_r = [birkach_r, stegaurach_r, strullendorf_west_r, strullendorf_nord_r]
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

# ------------------------------------------------------------------------------------
# Draw samples of L and x and pop if L < x
def sample_x_L(seed, number_of_samples, threshold, path, add=None):
    """
    Parameters
    ----------

    seed : integer
        Seed of random sampling.
    number_of_samples : integer
        Number of samples. Afterwards all sample combinations which don't meet the condition L > x are deleted!
    threshold : float
        Maximum value of L = aquifer length
    path : string
        Path to output directory.
    add : list of tuples [(L1, x1), (L2, x2), ...]
        Add your own values for L and x.
    """

    np.random.seed(seed)

    if rank == 0:
        # make output path
        if not os.path.exists(path):
            os.mkdir(path)

    Ls = []
    xs = []
    # epsilon = 5
    def get_samples(a, b, c, which="uni"):
        if which == "uni":
            samples = np.random.uniform(a, b, c)
        elif which == "lognorm":
            # a = mean
            # b = sigma
            # c = size
            samples = np.random.lognormal(a, b, c)
        else:
            print("Which mus be either 'uni' or 'log'.")
        return samples

    def erase_x_greater_L(Ls_sample, xs_sample):
        Ls = []
        xs = []
        for L, x in zip(Ls_sample, xs_sample):
            if (L > x) & (L < threshold):
                Ls.append(L)
                xs.append(x)
        return Ls, xs

    def normalize_by_mean(series):
        series_mean = np.mean(series)
        series_norm = [i / series_mean for i in series]
        return list(series_norm)

    # Ls_sample = get_samples(10, 40000, 200)
    # xs_sample = get_samples(10, 40000, 200)
    Ls_sample = get_samples(10, 4, number_of_samples, "lognorm")
    xs_sample = get_samples(10, 4, number_of_samples, "lognorm")
    Ls, xs = erase_x_greater_L(Ls_sample, xs_sample)
    Ls = Ls + [i[0] for i in add]
    xs = xs + [i[1] for i in add]
    Ls = np.around(Ls, 0)
    xs = np.around(xs, 0)
    L_and_x_input = [(L, x) for L, x in zip(Ls, xs)]
    # L_and_x_input.append([(L, x) for L, x in add])
    L_and_x_input.sort()

    # ------------------------------------------------------------------------------------
    # plot the samples
    if rank == 0:
        np.savetxt(path + "/" + comment + "L_and_x_input.txt", L_and_x_input)
        # Scatter plot
        plt.figure(figsize=(20, 10))
        plt.scatter(
            x=[i[0] for i in L_and_x_input], y=[i[1] for i in L_and_x_input], s=10
        )
        plt.yscale("log")
        plt.xscale("log")
        plt.ylabel("x")
        plt.xlabel("L")
        plt.ylim((10, 55000))
        plt.xlim((10, 55000))
        plt.title("Input parameter for L and x")
        plt.savefig(path + "/" + comment + "L_and_x_input_scatter.png", dpi=300)
        plt.close()
        # line plot
        plt.figure(figsize=(20, 10))
        L_plot, x_plot = plt.plot(L_and_x_input)
        plt.ylabel("[m]")
        plt.xlabel("count")
        plt.title("Input parameter for L and x")
        plt.legend([L_plot, x_plot], ["aquifer length", "position"], loc=1)
        plt.savefig(
            path + "/" + comment + "L_and_x_input_plot.png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()
        # histogram
        plt.figure(figsize=(20, 10))
        plt.hist(Ls, bins=20, alpha=0.5, label="aquifer length")
        plt.hist(xs, bins=20, alpha=0.5, label="position")
        plt.ylabel("count")
        plt.xlabel("value")
        plt.title("Histograms of L and x input parameter")
        plt.legend()
        plt.savefig(
            path + "/" + comment + "L_and_x_input_hist.png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()
    return L_and_x_input, Ls, xs


# create a list of tuples [(L1, x1), (L2, x2), ...] with own values for L and x which should be included in the sampling distribution
own_list = [(2528, 304), (2405, 2100), (968, 220),  (6032, 4100), (17102,12000), (4619,4500), (7910,7250), (8100, 5600)]

"""
taken from 20200708_appl_sa_main_new.py
# set gis nearest
xs = [304, 2100, 220, 4100]
Ls = [2528, 2405, 968, 6032]
# set gis farthest
xs = [12000, 4500, 7250, 5600]
Ls = [17102, 4619, 7910, 8100]
paths_h = [birkach_h, stegaurach_h, strullendorf_west_h, strullendorf_nord_h]
paths_r = [birkach_r, stegaurach_r, strullendorf_west_r, strullendorf_nord_r]
"""

# seet sampling parameters
# good seeds: 95321, 8465312, 97864334
seed = 132456
number_of_samples = 1500
threshold = 45000
path = os.path.join(output_path, "sampling")
# sample!
L_and_x_input, Ls, xs = sample_x_L(
    seed, number_of_samples, threshold, path, add=own_list
)

# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------
# print("# ------------------------------------------------------------------------------------")
# print("TEMP DELETE !!!!!! ONLY FOR STUDY 20200813 to get frequencies for each obs")

# L_and_x_input = (1000, 500)
# Ls = [1000]
# xs = [500]
# print("# ------------------------------------------------------------------------------------")
# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# set path to results incl file name of results
path_to_results_df = os.path.join(
    output_path, comment + "results_from_rank" + str(rank) + ".csv"
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
