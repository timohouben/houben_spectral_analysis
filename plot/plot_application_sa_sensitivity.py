# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------
# import modules
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import string
# set matplotlib style
CWD = os.getcwd()
CWD_p = os.path.dirname(CWD)
plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
# ------------------------------------------------------------------------------
"""
This module contains different functions to plot the results for the spectral 
analysis of real groundwater head data considering different L and x. 
The input file of these functions is a .csv with the following header:
index,name_h,name_r,L,x,S_out,T_out,tc_out,cov,time_step_size,cut_freq_higher,cut_freq_lower,rsquared,rsquared_log,rmse,rmse_log,combine_df,detrend_ts
"""


def plot_scatters_for_variables(
    path,
    filename,
    names_h,
    xy,
    variables,
    scales,
    v_names=None,
    label_name_h=None,
    own_list=None,
    ext=".pdf"
):
    """
    Create a figure with several scatter plots representing the samples and the
    resulting values for each variable.
    
    Parameters
    ----------
    path : string
        Path to the results.csv
    filename : string
        Filename of the results.csv
    names_h : list of strings
        Names of the observation points like they are named in the results.csv
    xy : tuple of 2 strings
        Pass the names of the columns to be used for x and y in the scatterplot.
        e.g.: xy = ("L", "x")
    variables : list of strings
        Variables to plot like they are named results.csv
    scales : list of strings
        Specify the scale of the colormap. Either "log" or "lin" for each variable.
    v_names : list of strings
        Display names of the variables.
    label_names : list of strings
        Display names of the observation points.
    own_list : list of tuples, optionally
        List of tuples with each tuple representing (L, x) which will be taken 
        to highlight certain points. Specify for each observation points two
        tuples with (Lmin, xmin),(Lmax, xmax) and concatinate them to one list,
        in the order of names_h, resulting in a list with 2*len(names_h) entries.
    ext : string, optionally
        File extension of the saved image.
  
    """
    # get the path to the file including the filename
    filepath = os.path.join(path, filename)
    # get pandas data frame
    results = pd.read_csv(filepath)
    # remove rows with erroneous values
    results = results.drop(results[results["cov"] == "[[inf inf]\n [inf inf]]"].index)
    # define columns and rows of subplots based in input lists
    nrow = len(names_h)
    ncol = len(variables)
    # BEGIN ---------------------
    # NEW definition of cmap
    # gist_rainbow
    cmap = mpl.cm.get_cmap("brg")
    cmap.set_over("cyan")
    cmap.set_under("orange")
    # END -------------------- 
    # markersize
    markersize = 20
    # define a counter for each plot
    count = 0
    # define a list for numbers for plots
    alpha_list = ["({})".format(i) for i in string.ascii_lowercase] + ["({})".format(i) for i in string.ascii_uppercase]

    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(ncol*5+1,nrow*4+1))
    for row, (name_h, label_name_h) in enumerate(zip(names_h, label_names_h)):
        for col, (variable, v_name, scale) in enumerate(
            zip(variables, v_names, scales)
        ):
            # get the right axes
            ax = axes[row, col]
            # get the right series
            data = results[results["name_h"] == name_h][variable]
            # --------------------
            # plot histogramms of data, comment to plot scatterplot correctly
            plt.figure(figsize=(12,7))
            if variable == "T_out" or variable == "S_out":
                data = np.log10(data)
                plt.xlabel("log " + v_name)
            else:
                plt.xlabel(v_name)
            if variable == "tc_out":
                print("Minimum tc of well: " + name_h + " is " + str(min(data)))
                data.where(data < 1000, np.nan, inplace=True)
            plt.ylabel("Count")
            plt.hist(data, bins=50)
            plt.savefig(os.path.join(path,"hists_" + name_h + "_" + variable + ext))
            plt.close()
            # --------------------
            # get x and y data
            x = results[results["name_h"] == name_h][xy[0]]
            y = results[results["name_h"] == name_h][xy[1]]
            if own_list is not None:
                # find index in series for L (assuming this value appears only once)
                idx_mask = (x == own_list[row * 2][0]) | (x == own_list[row * 2 + 1][0])
                idx = [markersize * 30 if i is True else markersize for i in idx_mask]
            else:
                idx = [markersize] * len(x)
            if scale == "log":
                # find min and max for all observation points
                minimum = results[variable].min()
                maximum = results[variable].max()
                sc = ax.scatter(
                    x,
                    y,
                    c=data,
                    s=idx,
                    marker="*",
                    # activate this line to scale cbar for only this obs point
                    #norm=colors.LogNorm(vmin=data.min(), vmax=data.max()),
                    norm=colors.LogNorm(vmin=minimum, vmax=maximum),
                    cmap=cmap,
                )
            if scale == "lin":
                # scale cbar for this observation point only
                if variable == "tc_out":
                    data_adj = data[data.between(data.quantile(0.15), data.quantile(0.85))]
                elif variable == "rsquared":
                    data_adj = data[data.between(data.quantile(0.15), data.quantile(1))]
                else:
                    data_adj = data[data.between(data.quantile(0.05), data.quantile(0.95))]
                # scale colorbar for all observation points
                # data_series = results[variable]
                # data_adj = data_series[data_series.between(data_series.quantile(0.15), data_series.quantile(0.85))]
                norm = mpl.colors.Normalize(vmin=data_adj.min(), vmax=data_adj.max())
                sc = ax.scatter(
                    x,
                    y,
                    c=data,
                    s=idx,
                    marker="*",
                    #vmax=data_adj.max(),
                    #vmin=data_adj.min(),
                    norm=norm,
                    cmap=cmap,
                )

            ax.annotate(alpha_list[count], (2*1e2, 5*1e4), fontsize=16)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(10e1, 10e4)
            ax.set_ylim(10e1, 10e4)
            ax.grid(b=True, which="minor", axis="x")
            if row == 0:
                if v_name is not None:
                    ax.set_title(v_name)
                else:
                    ax.set_title(variable)
            if row == nrow - 1:
                ax.set_xlabel("L [m]")
            else:
                ax.set_xticklabels([])
            if col == 0:
                if label_name_h is not None:
                    ax.set_ylabel(label_name_h + "\nx [m]")
                else:
                    ax.set_ylabel(name_h + "\nx [m]")
            else:
                ax.set_yticklabels([])
            fig.colorbar(sc, ax=axes[row, col], extend="max")
            count += 1
    plt.savefig(path + "scatters_variables_" + "_".join(variables) + ext)
    plt.show()

def plot_spectra_for_samples():
    """
    Plots a power spectrum for each observation point and the fittet analytical
    solution 
    """


if __name__ == "__main__":
    pass

    # plot every variable
    # --------------------------------------------------------------------------
    # the path to the file (upper: old, lower: new inlcusive edited L and x of Birkach)
    #path = "/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/output/combined_results/"
    path = "/Users/houben/phd/studies/application_spectral_analysis/main/20200729_spectral_analysis_sensitivity/20210708_edit_l_x_birkach/"
    # the file name (upper: old, lower: new inlcusive edited L and x of Birkach)
    #filename = "1_results_from_rank36_merge.csv"
    filename = "1_results_from_rank36_merge_inkl_edit_birkach.csv"
    # list of names of observations from column name_h
    names_h = ["Birkach", "Stegaurach", "Strullendorf_West", "Strullendorf_Nord"]
    label_names_h = ["Birkach", "Stegaurach", "Strullendorf West", "Strullendorf Nord"]
    # list of variables to be plotted from column header
    variables = [
        "S_out",
        "T_out",
        "tc_out",
        "rsquared",
        "rsquared_log",
        "rmse",
        "rmse_log",
    ]
    # list of names of variables used for title
    v_names = [
        "$S\ [-]$",
        "$T\ [m^2s^{-1}]$",
        "$t_c\ [days]$",
        "$R^2 [-]$",
        "$R^2\ log [-]$",
        "$RMSE$",
        "$RMSE\ log$",
    ]
    # define the colorbar scale for each variable
    scales = ["log", "log", "lin", "lin", "lin", "lin", "lin"]
    # provide list with manually added points
    own_list = [
        (2528, 304),
       # (17102, 12000), # old ones from Birkach
        (7509, 2997), # new ones form Birkach
        (2405, 2100),
        (4619, 4500),
        (968, 220),
        (7910, 7250),        
        (6032, 4100),
        (8100, 5600),
    ]
    # define x and y
    xy = ("L", "x")

    # plot_scatters_for_variables(
    #     path=path,
    #     filename=filename,
    #     names_h=names_h,
    #     label_name_h=label_names_h,
    #     variables=variables,
    #     v_names=v_names,
    #     scales=scales,
    #     own_list=own_list,
    #     xy=xy,
    # )
    
    # plot only S, T, tc and R2 for the paper
    # --------------------------------------------------------------------------
    
    variables = [
        "S_out",
        "T_out",
        "tc_out",
        "rsquared",
    ]

    plot_scatters_for_variables(
    path=path,
    filename=filename,
    names_h=names_h,
    label_name_h=label_names_h,
    variables=variables,
    v_names=v_names,
    scales=scales,
    own_list=own_list,
    xy=xy,
    )

    # extract some values for each observation well
    # --------------------------------------------------------------------------
    results = pd.read_csv(os.path.join(path, filename))
    birkach = results[results.name_h == "Birkach"]
    stegaurach = results[results.name_h == "Stegaurach"]
    strull_w = results[results.name_h == "Strullendorf_West"]
    strull_n = results[results.name_h == "Strullendorf_Nord"]
    
    # rows with max r2
    birkach[birkach.rsquared == birkach.rsquared.max()]
    birkach[birkach.rsquared_log == birkach.rsquared_log.max()]
    stegaurach[stegaurach.rsquared == stegaurach.rsquared.max()]
    strull_w[strull_w.rsquared == strull_w.rsquared.max()]
    strull_n[strull_n.rsquared == strull_n.rsquared.max()]

    # S for man samp
    birkach.S_out[birkach.L == own_list[0][0]]
    birkach.S_out[birkach.L == own_list[1][0]]
    stegaurach.S_out[stegaurach.L == own_list[2][0]]
    stegaurach.S_out[stegaurach.L == own_list[3][0]]
    strull_w.S_out[strull_w.L == own_list[4][0]]
    strull_w.S_out[strull_w.L == own_list[5][0]]
    strull_n.S_out[strull_n.L == own_list[6][0]]
    strull_n.S_out[strull_n.L == own_list[7][0]]
    # T for man samp
    birkach.T_out[birkach.L == own_list[0][0]]
    birkach.T_out[birkach.L == own_list[1][0]]
    stegaurach.T_out[stegaurach.L == own_list[2][0]]
    stegaurach.T_out[stegaurach.L == own_list[3][0]]
    strull_w.T_out[strull_w.L == own_list[4][0]]
    strull_w.T_out[strull_w.L == own_list[5][0]]
    strull_n.T_out[strull_n.L == own_list[6][0]]
    strull_n.T_out[strull_n.L == own_list[7][0]]
    # tc for man samp
    birkach.tc_out[birkach.L == own_list[0][0]]
    birkach.tc_out[birkach.L == own_list[1][0]]
    stegaurach.tc_out[stegaurach.L == own_list[2][0]]
    stegaurach.tc_out[stegaurach.L == own_list[3][0]]
    strull_w.tc_out[strull_w.L == own_list[4][0]]
    strull_w.tc_out[strull_w.L == own_list[5][0]]
    strull_n.tc_out[strull_n.L == own_list[6][0]]
    strull_n.tc_out[strull_n.L == own_list[7][0]]
    # rsquared for man samp
    birkach.rsquared[birkach.L == own_list[0][0]]
    birkach.rsquared[birkach.L == own_list[1][0]]
    stegaurach.rsquared[stegaurach.L == own_list[2][0]]
    stegaurach.rsquared[stegaurach.L == own_list[3][0]]
    strull_w.rsquared[strull_w.L == own_list[4][0]]
    strull_w.rsquared[strull_w.L == own_list[5][0]]
    strull_n.rsquared[strull_n.L == own_list[6][0]]
    strull_n.rsquared[strull_n.L == own_list[7][0]]