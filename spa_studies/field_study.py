# -*- coding: utf-8 -*
"""
Study to perform the spectral analysis on 4 groundwater wells in the Main
catchment. This study has been performed previously but with time series
which ended in 2009. This study has been extended to 2019 with new data from
mHM model run until 2019.

The following steps are perfomed:

1) The data is loaded and preprocessed
2) The spectra is calculated
3) The spectra is fitted with the analytical solution
4) The spectra are saved as images and txt
5) Results images are generated and saved. A resulting dataframe is saved
    with all parameters for each fit.

The spectral analysis will be performed with different parameters:

- detrend: yes/no
- cutted spectrum: 10e-7 ...
- different L and x (use specific)

The period with maximum overlap will be taken to calculate the spectra.
"""
# -----------------------------------------------------------------------------
# import modules
import sys
import numpy as np
import pandas as pd
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
from analysis.evaluate_results import calculate_rsquared
# path to save images
save_path = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/output"
# Load head
birkach_h = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/gwm/birkach_quadratic.txt"
stegaurach_h = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/gwm/stegaurach_quadratic.txt"
strullendorf_west_h = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/gwm/strull_west_quadratic.txt"
strullendorf_nord_h = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/gwm/strull_nord_quadratic.txt"
# load recharge
birkach_r = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/recharge/4416000.0_5522000.0_birkach_recharge_2019.txt"
stegaurach_r = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/recharge/4416000.0_5526000.0_stegaurach_recharge_2019.txt"
strullendorf_west_r = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/recharge/4424000.0_5526000.0_strullendorf_nord_recharge_2019.txt"
strullendorf_nord_r = "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/recharge/4424000.0_5522000.0_strullendorf_west_recharge_2019.txt"
names_h = ["Birkach", "Stegaurach", "Strullendorf_West", "Strullendorf_Nord"]
names_r = ["Birkach", "Stegaurach", "Strullendorf_West", "Strullendorf_Nord"]
# -----------------------------------------------------------------------------
# Set the follwing parameters manually
# set A alt
xs = [25800, 600, 7000, 1000]
Ls = [35000, 3000, 8000, 900]
# set B alt
xs = [2580, 6000, 70000, 10000]
Ls = [3500, 30000, 80000, 9000]
# set farthest
xs = [27500, 4700, 7000, 36000]
Ls = [37000, 4500, 7860, 38500]
# set nearest
xs = [500, 1700, 150, 5600]
Ls = [2800, 2000, 800, 8100]
# set gis nearest
xs = [304, 2100, 220, 4100]
Ls = [2528, 2405, 968, 6032]
# set gis farthest
xs = [12000, 4500, 7250, 5600]
Ls = [17102, 4619, 7910, 8100]
paths_h = [birkach_h, stegaurach_h, strullendorf_west_h, strullendorf_nord_h]
paths_r = [birkach_r, stegaurach_r, strullendorf_west_r, strullendorf_nord_r]
m = None
n = None
norm = False
convergence = 0.1
time_step_size = 86400
# cut higher frequencies than cut_freq_higher
cut_freq_higher = 1e-6
# cut lower frequencies than cut_freq_lower
cut_freq_lower = 0
combine_df = True
detrend_ts = True

# -----------------------------------------------------------------------------
for name_h, path_h, name_r, path_r, x, L in zip(names_h, paths_h, names_r, paths_r, xs, Ls):
    print("Currently fitting: " + str(name_h) + " ...")
    # Load the data
    head_df = pd.read_csv(
        path_h,
        sep=" ",
        header=None,
        names=["date", "head"],
    )
    head_df["date"] = pd.to_datetime(head_df["date"])

    # load recharge data
    recharge_df = pd.read_csv(
        path_r,
        sep=" ",
        header=None,
        names=["date", "recharge"],
    )
    recharge_df["date"] = pd.to_datetime(recharge_df["date"])


    if combine_df is True:
        # combine dataframes and remove rows with nans
        combined_df = pd.merge_ordered(recharge_df, head_df, how="inner")
        date_min = combined_df["date"].min()
        date_max = combined_df["date"].max()
        period = combined_df["date"].max() - combined_df["date"].min()
        print("Start/end/length of series where head measurements and recharge overlap: " + str(date_min) + "/" + str(date_max) + "/" + str(period))
        recharge_time_series = combined_df["recharge"].tolist()
        head_time_series = combined_df["head"].tolist()
    else:
        recharge_time_series = recharge_df["recharge"].tolist()
        head_time_series = head_df["head"].tolist()
        # modify the time series so that both have same length
        # assume: recharge is longest
        recharge_time_series = recharge_time_series[-len(head_time_series):]

    if detrend_ts is True:
        head_time_series = detrend(head_time_series)
        recharge_time_series = detrend(recharge_time_series)
    else:
        pass

    # convert mm/d to recharge along the aquifer in m2/s
    recharge_time_series = [i / 86400 / 1000 for i in recharge_time_series]

    '''
    ####################################################################
    # artificial data to test the script
    # A) S = 0.5, T = 0.008???, L = 1000, x = 200, white noise
    recharge_time_series = np.loadtxt("/Users/houben/phd/modelling/20190318_spectral_analysis_homogeneous/models_test/1001_24.1127_5.00e-01_8.00e-02/rfd_curve#1.txt")
    head_time_series = np.loadtxt("/Users/houben/phd/modelling/20190318_spectral_analysis_homogeneous/models_test/1001_24.1127_5.00e-01_8.00e-02/head_ogs_obs_00200_mean.txt")
    x = 200
    # B) S = 1.1e-5, T = 1.0 e-3, L = 1000, x = 200, white noise
    recharge_time_series = np.loadtxt("/Users/houben/phd/modelling/20190304_spectral_analysis_homogeneous/models/100_sample2_351_1.10e-05_1.00e-03/rfd_curve#1.txt")
    head_time_series = np.loadtxt("/Users/houben/phd/modelling/20190304_spectral_analysis_homogeneous/models/100_sample2_351_1.10e-05_1.00e-03/head_ogs_obs_00200_mean.txt")
    x = 200
    # C) 1076_border_50_stor_0.0001_rech_mHM, S = 1e-4, T2 = 3e-2
    #recharge_time_series = np.loadtxt("/Users/houben/phd/modelling/20190717_SA_hetero_block_2/1076_border_50_stor_0.0001_rech_mHM/rfd_curve#1_y_values.txt")
    #head_time_series = np.loadtxt("/Users/houben/phd/modelling/20190717_SA_hetero_block_2/1076_border_50_stor_0.0001_rech_mHM/head_ogs_obs_00500_mean.txt")
    #x = 500
    ############
    L = 1000
    time_step_size = 86400
    convergence = 0.01
    m = None
    n = None
    comment = "3_"
    norm = False
    # limits for the spectrum plot
    lims_head = [(1e-9,6e-6),(1e-6,1e12)]
    lims_base = [(1e-9,6e-6),(1e-6,3e1)]
    # cut the data
    # cut_value of frequency
    begin_index = 0
    end_index = 10000
    shift = 0
    #%matplotlib qt
    import matplotlib.pyplot as plt
    #plt.plot(recharge_time_series)
    #plt.show()
    recharge_time_series = recharge_time_series[begin_index:end_index]
    head_time_series = head_time_series[begin_index+shift:end_index+shift]
    ####################################################################
    '''

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


    # fit the power spectrum with the analytical solution
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
        popt, pcov = [np.nan, np.nan], [[np.nan, np.nan],[np.nan, np.nan]]
        print("popt and pcov have been set to np.nan")
    except ValueError:
        print("either ydata or xdata contain NaNs, or if incompatible options are used")
        popt, pcov = [np.nan, np.nan], [[np.nan, np.nan],[np.nan, np.nan]]
    except OptimizeWarning:
        print("Covariance of the parameters could not be estimated.")
        #popt, pcov = [np.nan, np.nan], [[np.nan, np.nan],[np.nan, np.nan]]

    Sy = abs(popt[0])
    T = abs(popt[1])

    # calculate the fitted power spectra
    Shh_fitted = shh_analytical(
        (frequency_input, Sww),
        Sy,
        T,
        x,
        L,
        m=n,
        n=m,
        norm=norm,
    )


    # get the characteristic time
    tc = calc_tc(L,Sy,T,which="dupuit")
    # define a (discharge constant)
    #a = np.pi ** 2 * T / (4 * L ** 2)
    # define tc (characteristic time scale)
    #tc = Sy / a

    data = np.vstack((Shh, Shh_fitted))

    labels = [
        "Shh numerical",
        "Shh fitted"
    ]
    linestyle = ["-", "-"]
    marker = ["", "d"]
    figtxt ="Derived Parameter:    S = %1.3e, T = %1.3e [m2/s], tc = %1.3e [d]\nInput Parameter:        L = %0.0f, x = %0.0f" % (
        Sy,
        T,
        tc,
        L,
        x
    )
    rsquared = calculate_rsquared(data[0], data[1])
    print("Currently plotting: " + str(name_h) + " ...")
    # plot only spectrum of Shh
    plot_spectrum(
        [Shh],
        frequency_output,
        heading="Shh - Head Power Spectrum " + name_h,
        labels=["Shh observed"],
        path=save_path,
        linestyle=[""],
        marker=["."],
        lims=[(2e-9,7e-6),(1e-5,1e7)],
        name="Shh_" + name_h,
        color="blue"
    )

    # plot only spectrum of Sww
    plot_spectrum(
        [Sww],
        frequency_input,
        heading="Sww - Recharge Power Spectrum  " + name_r,
        labels=["Sww"],
        path=save_path,
        linestyle=[""],
        marker=["."],
        lims=[(2e-9,7e-6),(1e-20,1e-9)],
        name="Sww_" + name_r
    )

    # plot Shh and the fitted spectrum
    plot_spectrum(
        data,
        frequency_input,
        heading="Shh - Head Power Spectrum " + name_h,
        labels=["Shh observed", "Shh fitted"],
        path=save_path,
        linestyle=["-", " "],
        marker=["", "."],
        figtxt=figtxt,
        lims=[(2e-9,7e-6),(1e-7,1e8)],
        name="Shh_fitted_" + name_h + "_R_" + name_r + "_L_" + str(L) + "_x_" + str(x),
        r2=rsquared
    )

    # save the spectra and frequency as txt
    np.savetxt(save_path + "/" + name_h + ".txt", Shh)
    np.savetxt(save_path + "/" + name_r + ".txt", Sww)
    np.savetxt(save_path + "/" + "Shh_fitted_" + name_h + "_R_" + name_r + "_L_" + str(L) + "_x_" + str(x) + ".txt", Shh_fitted)
    np.savetxt(save_path + "/" + "frequency" + ".txt", frequency_output)
