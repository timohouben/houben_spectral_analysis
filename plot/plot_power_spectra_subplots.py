# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import string

# set matplotlib style
CWD = os.getcwd()
CWD_p = os.path.dirname(CWD)
plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])


def plot_head_and_recharge_power_spectrum_subplots(
    savepath, titles, labels, paths_y, paths_x, rows, cols, markers, linestyles, colors
):
    """
    Plot multiple power spectra in subplots.

    Parameters
    ----------
    savepath : string
        Path to save the plot.
    titles : list of strings
        Titles of the subplots.
    labels : list of strings
        Labels of the plots.
    paths_y : list of strings
        Paths of y data.
    paths_y : list of strings
        Paths of x data
    rows : list of integers
        Row of each plot for subplots.
    cols : list of integers
        Column of each plot for subplots
    markers : list of strings
        Markerstyle for plots.
    linestyles : list of strings
        Linestyles for plots.
    colors : list of strings
        Colors for plots
    
    Yields
    ------

    A .pdf plot in the savepath directory.

    """
    # define a list for numbers for plots
    alpha_list = ["({})".format(i) for i in string.ascii_lowercase] + [
        "({})".format(i) for i in string.ascii_uppercase
    ]
    # define frequencies for year, month, week
    year = 1 / 365 / 86400
    month = 1 / 30 / 86400
    week = 1 / 7 / 86400
    # get length of rows/cols
    nrows = len(np.unique(rows))
    ncols = len(np.unique(cols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharey=True,
        sharex=True,
        figsize=(ncols * 8 + 1, nrows * 4 + 1),
    )
    # iterate over input lists
    for title, label, path_y, path_x, row, col, marker, linestyle, color in zip(
        titles, labels, paths_y, paths_x, rows, cols, markers, linestyles, colors
    ):
        ax = axes[row, col]
        # load data
        x = np.loadtxt(path_x)
        y = np.loadtxt(path_y)
        # Status
        print(title)
        # plot
        ax.loglog(
            x, y, label=label, marker=marker, ls=linestyle, color=color, markersize=4
        )
        # make vertical lines for year, month, week
        for cycle, label in zip([year, month, week], ["Year", "Month", "Week"]):
            ax.vlines(
                cycle, ymin=0, ymax=1e8, linestyles="dotted", alpha=0.5, colors="grey"
            )
            if row == ncols - 1:
                ax.annotate(
                    label, (cycle * 1.1, 1e-5), alpha=0.5, fontsize="14", color="grey"
                )
        ax.set_title(title)
        ax.legend()
        if col == 0:
            ax.set_ylabel("Power Spectrum")
        if row == ncols - 1:
            ax.set_xlabel("Frequency [Hz]")

    # make an annotation for numbering of subplots
    for i, ax in enumerate(np.ravel(axes)):
        ax.annotate(alpha_list[i], (1e-9, 1e7), fontsize=16)
    plt.subplots_adjust(wspace=0.1)
    plt.savefig(
        os.path.join(
            savepath, "_".join([i.replace(" ", "") for i in titles[0::2]]) + ".pdf"
        )
    )


if __name__ == "__main__":

    # Configuration for the paper spectral_analysis_1
    # define mandatory input lists
    # ------------------------------------------------------------------------------
    savepath = "/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/output"
    # list with names
    titles = [
        "Birkach",
        "Birkach",
        "Stegaurach",
        "Stegaurach",
        "Strullendorf West",
        "Strullendorf West",
        "Strullendorf Nord",
        "Strullendorf Nord",
    ]

    labels = ["$S_{hh}$", "$S_{hh}\ fit$"] * 4

    # list with paths to y data
    paths_y = [
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/birkach/Shh_Birkach.txt",
        # best r2 from man samp
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/birkach/Shh_fitted_Birkach_R_Birkach_L_2528.0_x_304.0.txt",
        
        # best r2
        #"/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/birkach/Shh_fitted_Birkach_R_Birkach_L_8194.0_x_0.0.txt",
        
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/stegaurach/Shh_Stegaurach.txt",
        # best r2 from man samp
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/stegaurach/Shh_fitted_Stegaurach_R_Stegaurach_L_2405.0_x_2100.0.txt",
        # best r2
        #"/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/stegaurach/Shh_fitted_Stegaurach_R_Stegaurach_L_13297.0_x_0.0.txt",
        
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_w/Shh_Strullendorf_West.txt",
        # best r2 from man samp
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_w/Shh_fitted_Strullendorf_West_R_Strullendorf_West_L_7910.0_x_7250.0.txt",
        # best r2
        #"/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_w/Shh_fitted_Strullendorf_West_R_Strullendorf_West_L_6574.0_x_6415.0.txt",
        
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_n/Shh_Strullendorf_Nord.txt",
        # best r2 from man samp
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_n/Shh_fitted_Strullendorf_Nord_R_Strullendorf_Nord_L_6032.0_x_4100.0.txt",
        # best r2
        #"/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_n/Shh_fitted_Strullendorf_Nord_R_Strullendorf_Nord_L_8194.0_x_0.0.txt",
    ]

    # list with paths to x data
    paths_x = [
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/birkach/Frequency_input_Birkach.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/birkach/Frequency_input_Birkach.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/stegaurach/Frequency_output_Stegaurach.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/stegaurach/Frequency_output_Stegaurach.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_w/Frequency_output_Strullendorf_West.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_w/Frequency_output_Strullendorf_West.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_n/Frequency_input_Strullendorf_Nord.txt",
        "/Users/houben/Desktop/eve_data/20200729_spectral_analysis_sensitivity/output/strullendorf_n/Frequency_input_Strullendorf_Nord.txt",
    ]

    # list with row values
    rows = [0, 0, 0, 0, 1, 1, 1, 1]

    # list with column values
    cols = [0, 0, 1, 1, 0, 0, 1, 1]

    markers = ["", ".", "", ".", "", ".", "", "."]

    linestyles = ["-", " ", "-", " ", "-", " ", "-", " "]

    colors = [
        "#E24A33",
        "#FBC15E",
        "#348ABD",
        "#FBC15E",
        "#988ED5",
        "#FBC15E",
        "#777777",
        "#FBC15E",
    ]

    plot_head_and_recharge_power_spectrum_subplots(
        savepath,
        titles,
        labels,
        paths_y,
        paths_x,
        rows,
        cols,
        markers,
        linestyles,
        colors,
    )
    # ------------------------------------------------------------------------------
