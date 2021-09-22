# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division
import sys
sys.path.append("/Users/houben/phd/python/scripts/spectral_analysis")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# set some parameters for the analysis manually
dpi = 300
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


def plot_eff_thickness_vs_connection_length_and_thickness(results, path, obs_locs, stors):
    """
    Plot the effectie aquifer thickness which has been calculated from the
    derived t_c, the real input hydraulic conductivity and the real input
    spec. storage.

    e.g. for analysis 20191115_ogs_homogeneous_thickness

    Parameter
    ---------
    results : Pandas Dataframe
        Results dataframe from spectral analysis.
    path : string
        Path which will be used to store the resulting images.
    obs_locs : list
        List of observation points numbering
    stors : list
        List of storativities.

    Yields
    ------

    An image in the path directory.
    """

    import matplotlib.pyplot as plt
    import pandas as pd
    from cycler import cycler

    colors = [
        "#FF0000",
        "#FF3A00",
        "#FF7800",
        "#FF7800",
        "#FFFB00",
        "#B9FF00",
        "#59FF00",
        "#00FF9B",
        "#00FFE4",
        "#00ECFF",
        "#00D1FF",
        "#0097FF",
        "#0055FF",
        "#0017FF",
        "#4900FF",
        "#8F00FF",
        "#FB00FF",
        "#FF009B",
        "#FF0049",
        "#FF001B",
        "#FF0000",
    ]
    plt.rc("axes", prop_cycle=(cycler("color", colors)))
    results = results[results.S_in / results.aquifer_thickness == 0.0001]
    results = results[~results.effective_thickness.isnull()]

    plt.figure(figsize=(16, 9))
    for thickness in sorted(results.aquifer_thickness.unique(), reverse=False):
        # plot the absolut effective thickness
        plt.plot(results.c_length_rel[results.aquifer_thickness == thickness], results.effective_thickness[results.aquifer_thickness == thickness], label=thickness, ls="", marker="*")
        plt.ylabel("absolut effecive thickness")
        plt.ylim(0,1100)
        # plot the relative effective thickness
        #plt.plot(results.c_length_rel[results.aquifer_thickness == thickness], results.effective_thickness[results.aquifer_thickness == thickness] / results.aquifer_thickness[results.aquifer_thickness == thickness], label=thickness, ls="", marker="*")
        #plt.ylabel("relative effecive thickness")
        #plt.ylim(0,1.1)

    plt.legend()
    plt.title("effective aquifer thickness derived from baseflow spectral analysis")
    plt.xlabel("connection length / aquifer thickness")
    plt.savefig(path + "/abs_eff_thickness_vs_connection_length_by_aquifer_thickness.png")
    plt.show()


def plot_k_vs_connection_length_and_thickness(results, path, obs_locs, stors):
    """
    Plot the derived kf value from the SA vs the ration between the connection
    ength of the aquifer and the river. Each graph represents a thickness.

    Parameter
    ---------
    results : Pandas Dataframe
        Results dataframe from spectral analysis.
    path : string
        Path which will be used to store the resulting images.
    obs_locs : list
        List of observation points numbering
    stors : list
        List of storativities.

    Yields
    ------

    An image in the path directory.
    """

    import matplotlib.pyplot as plt
    import pandas as pd
    from cycler import cycler

    colors = [
        "#FF0000",
        "#FF3A00",
        "#FF7800",
        "#FF7800",
        "#FFFB00",
        "#B9FF00",
        "#59FF00",
        "#00FF9B",
        "#00FFE4",
        "#00ECFF",
        "#00D1FF",
        "#0097FF",
        "#0055FF",
        "#0017FF",
        "#4900FF",
        "#8F00FF",
        "#FB00FF",
        "#FF009B",
        "#FF0049",
        "#FF001B",
        "#FF0000",
    ]
    plt.rc("axes", prop_cycle=(cycler("color", colors)))

    results["k_out"] = results["T_out"] / results["aquifer_thickness"]
    results["Ss_in"] = results["S_in"] / results["aquifer_thickness"]
    k_in = (results.T_in / results.aquifer_thickness).unique()
    for stor in stors:
        for obs_loc in obs_locs:
            plt.figure(figsize=(16, 9))
            for aquifer_thickness in sorted(results.aquifer_thickness.unique()):
                # print(aquifer_thickness)
                # print("länge: " + str(len(results["k_out"][(results["aquifer_thickness"] == aquifer_thickness) & (results["obs_loc"] == obs_locs) & (results["Ss_in"] == 0.0001)])))
                # print("name: " + results["name"][(results["aquifer_thickness"] == aquifer_thickness) & (results["obs_loc"] == obs_locs) & (results["Ss_in"] == 0.0001)])
                if (
                    len(
                        results["k_out"][
                            (results["aquifer_thickness"] == aquifer_thickness)
                            & (results["obs_loc"] == obs_loc)
                            & (results["Ss_in"] == stor)
                        ]
                    )
                    == 10
                ):
                    pass
                    plt.semilogy(
                        results["c_length_rel"][
                            (results["aquifer_thickness"] == aquifer_thickness)
                            & (results["obs_loc"] == obs_loc)
                            & (results["Ss_in"] == stor)
                        ],
                        results["k_out"][
                            (results["aquifer_thickness"] == aquifer_thickness)
                            & (results["obs_loc"] == obs_loc)
                            & (results["Ss_in"] == stor)
                        ],
                        ls=" ",
                        marker="*",
                        label=aquifer_thickness,
                    )
            plt.title("Derived K for observation point " + str(obs_loc) + " and different aquifer thicknesses. Specific Storage: " + str(stor))
            plt.xlabel("Connection length / aquifer thickness")
            plt.ylabel("K [m/s]")
            plt.hlines(y=k_in, xmin=0, xmax=1)
            plt.ylim(1e-6, 1e-4)
            plt.legend(title="thickness")
            # plt.show()
            plt.savefig(
                path + "/" + str(obs_loc) + "_stor_" + str(stor) + "_k_out_vs_lc.png",
                dpi=300,
            )
            plt.close()


def plot_baseflow_sa_vs_boundary_block(results, path_to_results):
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    results_temp = results[
        results["recharge"]
        == "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"
    ]
    results_temp = results_temp[results_temp["S_in"] == 0.003]
    results_temp = results_temp[results_temp["obs_loc"] == 1000]
    results_temp = results_temp.sort_values(by="name")
    boundary = [
        50,
        100,
        150,
        200,
        250,
        300,
        350,
        400,
        450,
        500,
        550,
        600,
        650,
        700,
        750,
        800,
        850,
        920,
        930,
        940,
        950,
        960,
        970,
        980,
        990,
    ]
    # plt.semilogy(boundary, results_temp["D"])
    plt.semilogy(boundary, [float(i[1:-1]) for i in results_temp["D_cov"]])
    plt.grid(which="both")
    plt.ylabel("covariance of optimization")
    # plt.ylabel("Diffusivity [m^2/s]")
    plt.xlabel("boundary")
    plt.show()


def plot_parameter_vs_location_block(
    results,
    path_to_results,
    borders,
    S_in,
    recharge_list,
    threshold=0.1,
    comment="",
    saveimg=True,
    facecolor="#C0C0C0"
):
    """
    Plot the derived parameter (T or S) along the different observation points
    and for a selection of positions of the boarder between a high and a low
    conductive zone.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results from spectral analysis.
    path_to_results : string
        Path where to store the images.
    borders : list of integers
        List of values for the boarder between two zones.
    S_in : float
        For which S input.
    recharge_list : list of strings
        Names of the recharges.
    threshold : float
        Threshold for variance from curve_fit. All other entries will be dropped from data frame.
    comment : string
        Give a commment.
    saveimg : bool, default: True
        True: Saves a .png image like it has been modified in the interactive backend. Will be saved after you have closed the window.
    facecolor : string, Default: "#C0C0C0"
        Facecolcor of the plot in HTML color format or as name.

    Yields
    ------

    A/multiple plot(s) in path_to_results directory.
    """

    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from analysis.processing import identify_numbers_from_string

    # set matplotlib style
    CWD = os.getcwd()
    CWD_p = os.path.dirname(CWD)
    plt.style.use(os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle"))

    # convert strings to numbers
    results["cov_numbers"] = results["cov"].apply(identify_numbers_from_string)
    results["sigma_S"] = results["cov_numbers"].apply(
        lambda x: float(x[0]) if x != [] else np.nan
    )
    results["sigma_T"] = results["cov_numbers"].apply(
        lambda x: float(x[3]) if x != [] else np.nan
    )
    # remove all lines where the covariance is higher than threshold
    results = results[results["sigma_T"] > -threshold]
    results = results[results["sigma_T"] < threshold]

    fig, axs = plt.subplots(len(borders), 1, sharex=True, figsize=(20, len(borders)+2))
    for i, border in enumerate(borders):
        results_temp = results
        border_str = "border_" + str(border) + "_"
        # get a new column with border values
        results_temp = results_temp[results_temp.name.str.contains(border_str)]
        # only values with specific S_in
        results_temp = results_temp[results_temp["S_in"] == S_in]
        # plot the T in
        T_in = [
            float(results_temp.T_in_1.unique())
            if obs < border
            else float(results_temp.T_in_2.unique())
            for obs in np.arange(0, 1010, 10)
        ]
        axs[i].semilogy(
            np.arange(0, border, 10),
            T_in[: len(np.arange(0, border, 10))],
            linestyle="--",
            linewidth=8,
            marker="",
            label="high conductive part",
            color="#ff7f0e",

        )
        axs[i].semilogy(
            np.arange(border, 1010, 10),
            T_in[len(np.arange(0, border, 10)) :],
            linestyle="--",
            linewidth=8,
            marker="",
            label="low conductive part",
            color="#1f77b4",
        )
        axs[i].grid(which="both", color="grey", linestyle=":")
        # for both recharges in one plot
        colors = [
            "#CC3810",
            "#23E1D0",
            ]
        markers = [
            "^",
            "*",
            ]
        markersizes = [
            "10",
            "5"
            ]
        linewidths = [
            "5",
            "1.5",
            ]
        labels = [
            "white noise",
            "mHM"
        ]
        for recharge, color, marker, markersize, linewidth, label in zip(recharge_list, colors, markers, markersizes, linewidths, labels):
            results_temp_r = results_temp[results_temp["recharge"] == recharge]
            print("Plotting " + recharge)
            axs[i].plot(
                results_temp_r.obs_loc,
                results_temp_r.T_out,
                label="derived Transmissivity: " + label,
                linewidth=linewidth,
                marker=marker,
                markersize=markersize,
                color=color,
            )
        # Remove horizontal space between axes
        fig.subplots_adjust(hspace=0.2)
        # axs[i].set_title("Border at " + str(border) + " m")
        axs[i].annotate("Boundary \n" + str(border), (border + 10, 5e-4), fontsize=20)
        axs[i].axvline(x=border, color="black")
        #axs[i].spines["right"].set_visible(False)
        #axs[i].spines["top"].set_visible(False)
        #axs[i].spines["bottom"].set_visible(False)
        #axs[i].spines["left"].set_visible(True)
        # axs[i].set_ylim(ymin=np.min(T_in + results_temp.T_out.tolist())*0.5, ymax=np.max(T_in + results_temp.T_out.tolist())*4)
        axs[i].set_ylim(bottom=7e-5, top=5e-1)
        axs[i].set_facecolor(facecolor)
        # set second y achsis
        # axs_twin = axs[i].twinx()
        # axs_twin.bar(results_temp.obs_loc, results_temp.sigma_T)
        # axs_twin.errorbar(results_temp.obs_loc, results_temp.T_out, results_temp.sigma_T, label="Border at " + str(border) + " m", marker="+")
    fig.text(0.04, 0.5, r"Transmissivity $[m^2s^{-1}]$", va="center", rotation="vertical", fontsize=22)
    plt.xlabel("location [m]")
    fig.suptitle("Derived Transmissivity vs input Transmissivity\n" + comment)
    plt.legend(loc="lower left")
    axs[len(borders) - 1].legend(
        loc="upper center", bbox_to_anchor=(0.47, -0.75), fancybox=True, shadow=False, ncol=5
    )
    if saveimg == True:
        plt.savefig(path_to_results + "/" + comment + "T_vs_location.pdf", dpi=dpi)
    #plt.show()
    return fig

def plot_error_vs_tc(
    results, path_to_results, lims=None, locs=None, comment="", abs=True, yaxis="log"
):
    """
    Plot errors of input and output parameters (S, T, tc) vs input tc.
    This results in 3 (S,T,tc) plots with multiple graphs according to the
    different observation points.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results generated by multi_psd.py.
    path_to_results : string
        Path where to store the images
    abs : bool
        Take the mean over absolute errors or not.
    locs : list of integers
        default: [200, 400, 600, 800, 990]
    lims : list of two tuples with limits for x and y. [(x_min, x_max), (y_min, y_max)]
        Default: No limits
    yaxis : string, Default: "log"
        "log": loglog plot
        "lin": semilogx
    """

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    from cycler import cycler

    try:
        if locs == None:
            locs = [200, 400, 600, 800, 990]
    except ValueError:
        pass

    plt.style.use("ggplot")
    # ['#f1eef6', '#d0d1e6', '#a6bddb', '#74a9cf','#3690c0','#0570b0','#034e7b']
    # plt.rc(
    #     "axes",
    #     prop_cycle=(
    #         cycler(
    #             "color",
    #             [
    #                 #    '#edf8b1',
    #                 "#c7e9b4",
    #                 #    '#7fcdbb',
    #                 "#41b6c4",
    #                 "#1d91c0",
    #                 "#225ea8",
    #                 "#253494",
    #             ],
    #         )
    #     ),
    # )
    font = {"family": "DejaVu Serif", "weight": "normal", "size": 14}
    plt.rc("font", **font)
    if not os.path.exists(path_to_results + "/error_vs_tc"):
        os.mkdir(path_to_results + "/error_vs_tc")
    for error in ["err_S", "err_T", "err_tc"]:
        tc_agg = results["tc_in"].apply(lambda x: np.around(x, 2)).unique()
        tc_agg.sort()
        for loc in locs:  # results.obs_loc.unique():
            results_loc = results[error][results["obs_loc"] == loc]
            err_vs_tc_at_loc = []
            for tc in tc_agg:
                print(
                    "...currently grouping values for tc = " + str(tc),
                    "location: " + str(loc),
                    "Error: " + error,
                )
                # append error to list for specific tc
                if abs == True:
                    err_vs_tc_at_loc.append(
                        np.mean(
                            np.abs(
                                results_loc[
                                    results["tc_in"].apply(lambda x: np.around(x, 2))
                                    == tc
                                ]
                            )
                        )
                    )
                if abs == False:
                    err_vs_tc_at_loc.append(
                        np.mean(
                            results_loc[
                                results["tc_in"].apply(lambda x: np.around(x, 2)) == tc
                            ]
                        )
                    )
            if yaxis == "log":
                plt.loglog(tc_agg, err_vs_tc_at_loc, label=str(loc) + " m")
            elif yaxis == "lin":
                plt.semilogx(tc_agg, err_vs_tc_at_loc, label=str(loc) + " m")
            else:
                print("Please set argument 'yaxis' either to 'lin' or 'log'.")
        if lims != None:
            plt.ylim(lims[1])
            plt.xlim(lims[0])
        if error == "err_S":
            plt.ylabel(r"Error S [%]")
        elif error == "err_T":
            plt.ylabel(r"Error T [%]")
        else:
            plt.ylabel(r"$Error\ t_c\ [%]$")
        plt.xlabel(r"$t_c\ [days]$")
        plt.legend()
        plt.savefig(
            path_to_results
            + "/error_vs_tc/"
            + comment
            + "_"
            + error
            + "_vs_tc_abs_"
            + str(abs)
            + ".eps",
            dpi=dpi,
            bbox_inches="tight",
        )
        plt.close()


def plot_errors_vs_loc_aggregate(
    results, path_to_results, error, aggregate, bins, abs, comment=""
):
    """
    Plot errors of input and output parameters (S, T, tc) vs the observation
    location in the aquifer (2D transect) and aggregated for specific ranges of tc.
    This results in three plots over all model runs.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results generated by multi_psd.py.
    path_to_results : string
        Path where to store the images
    error : string
        Select the column from data frame to plot along locations.
    aggregate : string
        Must be an existing column in dataframe. Lines in plot will be aggregated based on the chosen column.
    bins : list of floats
        Lines in plot will be aggregated according to the values in the list and the chosen column.
    abs : bool
        True: Neglect the sign of errors an avegage only positive deviations.
    comment : string
        Give a meaningful comment for the file name. Otherwise it will override existing files.


    Yields
    ------

    see docstring description
    """
    # print(results)
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    print("Your have provided the following bins: ", bins)
    # get list of locations from df
    obs_locs = results["obs_loc"].unique()
    obs_locs = [int(i) for i in obs_locs]
    # make arrays with errors for obs_locs
    err = np.zeros((len(bins), len(obs_locs)))
    for i, loc in enumerate(obs_locs):
        for j, bin in enumerate(bins):
            if abs == True:
                temp = results[results["obs_loc"] == loc][results[aggregate] <= bin][
                    error
                ]
                temp.tolist()
                temp = [np.abs(i) for i in temp]
                err[j, i] = np.mean(temp)
            else:
                err[j, i] = np.mean(
                    results[results["obs_loc"] == loc][results[aggregate] <= bin][error]
                )

    for i in np.arange(len(bins)):
        plt.plot(obs_locs, err[i, :], label="< " + "{0:0.3e}".format(bins[i]))
    plt.legend()
    plt.title(("Plot: " + error + ", " + "Aggregation: " + aggregate))
    plt.ylabel(error)
    plt.xlabel("location [m]")
    plt.savefig(
        path_to_results + "/" + error + "-" + aggregate, bbox_inches="tight", dpi=dpi
    )
    plt.close()


def plot_errors_vs_loc_hetero(obs, error_list, legend, ylabel, path, comment):
    """
    Plot the error vs location in the aquifer.

    Parameters
    ----------

    obs : list
        List of x value of observations points
    error_list : list of lists for different errors
        List of errors for each observation point.
    legend : list of strings for different errors
        String for ylabel
    ylabel : string
        ylabel
    path : string
        Path where to srote the images.
    comment : string
        Give a comment to be stored in the filename.

    Yields
    ------

    Saves a plot in the path directory.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import os.path

    for i, error in enumerate(error_list):
        plt.plot(obs, error, label=legend[i])
    plt.hlines(0, 0, np.max(obs), colors="k", linestyles="dashed")
    plt.legend()
    plt.ylabel(ylabel)
    plt.title("Error vs location: " + os.path.basename(path))
    plt.savefig(path + "/" + comment + "Error_vs_loc", dpi=dpi, bbox_inches="tight")


def plot_errors_vs_loc_homo(results, path_to_results, comment=""):
    """
    Plot errors of input and output parameters (S, T, tc) vs the observation
    location in the aquifer (2D transect). This results in three plots per OGS
    model run.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results generated by multi_psd.py.
    comment : string
        Give a meaningful comment for the file name. Otherwise it will override existing files.

    Yields
    ------

    see docstring description
    """

    import matplotlib.pyplot as plt

    def plot(x, y, path, label, title):
        plt.plot(x, y, label=label)
        plt.xlabel("Location of observation point [m]")
        plt.ylabel("Error in %")
        plt.title(title)
        plt.legend()
        plt.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close()

    for project_folder in results.name.unique():
        err_S = results[results["name"] == project_folder]["err_S"]
        err_T = results[results["name"] == project_folder]["err_T"]
        err_tc = results[results["name"] == project_folder]["err_tc"]
        obs_loc = results[results["name"] == project_folder]["obs_loc"]
        plot(
            obs_loc,
            err_S,
            path_to_results + "/" + comment + project_folder + "_err_S.png",
            label=project_folder,
            title="Relative error in storativity",
        )
        plot(
            obs_loc,
            err_T,
            path_to_results + "/" + comment + project_folder + "_err_T.png",
            label=project_folder,
            title="Relative error in transmissivity",
        )
        plot(
            obs_loc,
            err_tc,
            path_to_results + "/" + comment + project_folder + "_err_tc.png",
            label=project_folder,
            title="Relative error in characteristic time",
        )


def plot_heatmap(results, path_to_results, abs=True, comment=""):
    """
    Plot errors of input and output parameters (S, T, tc) vs the parameter
    range as heatmap. This results in three plots per location.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results generated by multi_psd.py.
    path_to_results : string
        Where to store the heatmaps.
    abs : bool
        Absolute error in %.
    comment : string
        Give a meaningful comment for the file name. Otherwise it will override existing files.

    Yields
    ------

    see docstring description
    """

    import seaborn as sns
    import os

    # extract input values for achsis limits
    # achsisticks_x = [1,5,10,15,20,25,30,35,40,45]
    # achsislabel_x = ["%1.2e" % i for i in achsisticks_x]
    # achsisticks_y = [1,5,10,15,20,25,30,35,40,45]
    # achsislabel_y = ["%1.2e" % i for i in achsisticks_y]

    def plot(pivotted, error):
        import numpy as np
        from matplotlib.colors import LogNorm
        import math

        # set axismin and axis max based on input space (hard coded, BAD SOLUTION)
        # achsismin_y, achsismax_y = 1e-6, 1e-1
        # achsismin_x, achsismax_x = 1e-5, 1
        # achsisticks_x = [math.pow(10, i) for i in range(math.floor(math.log10(achsismin_x)), 1+math.ceil(math.log10(achsismax_x)))]
        # achsisticks_y = [math.pow(10, i) for i in range(math.floor(math.log10(achsismin_y)), 1+math.ceil(math.log10(achsismax_y)))]
        barmin, barmax = 1, 1000
        cbar_ticks = [1, 10, 100, 1000]
        log_norm = LogNorm(vmin=barmin, vmax=barmax)
        ax = sns.heatmap(
            pivotted,
            cmap="PuBu",
            cbar_kws={"ticks": cbar_ticks},
            vmax=barmax,
            vmin=barmin,
            norm=log_norm,
        )  # , yticklabels=achsislabel_y, xticklabels=achsislabel_x)
        # cmap="Spectral_r",
        # ax.set_yticks(achsisticks_y)
        # ax.set_xticks(achsisticks_x)
        ax.invert_yaxis()
        # import matplotlib.ticker as ticker
        # tick_locator = ticker.MaxNLocator(12)
        # ax.xaxis.set__locator(tick_locator)
        # ax.yaxis.set_major_locator(tick_locator)

        fig = ax.get_figure()
        # fig.set_size_inches(5, 5)
        if not os.path.exists(path_to_results + "/" + comment + "heatmap"):
            os.mkdir(path_to_results + "/" + comment + "heatmap")

        fig.savefig(
            path_to_results
            + "/"
            + comment
            + "heatmap"
            + "/"
            + str(obs_loc)
            + "_"
            + error,
            dpi=dpi,
            bbox_inches="tight",
        )
        fig.clf(fig)

    for obs_loc in results["obs_loc"].unique():
        # extract only rows with obs_loc==obs_loc
        df_obs_loc = results[results.obs_loc == obs_loc]
        # extract columns for plotting
        for error in ["err_S", "err_T", "err_tc"]:
            # absolute erros, NOT A GOOD SOLUTION
            if abs == True:
                results[error] = results[error].apply(
                    lambda x: x * (-1) if x < 0 else x
                )
            df_obs_loc_cut = df_obs_loc[["S_in", "T_in", error]]
            # pivot this table
            pivot_df_obs_loc_cut = df_obs_loc_cut.pivot("S_in", "T_in", error)
            # plot heatmap
            plot(pivot_df_obs_loc_cut, error)


def plot_parameter_vs_location(
    path_to_results, parameter, location, y_label, error=None
):
    """
    Plots a list/array of parameters along location and saves it in path_to_results

    Parameters
    ----------

    path_to_results : string
        Path to results.
    parameter : 1D List, array
        Parameter to plot
    location : 1D list, array
        Locations
    y_label : string
        Give the y-achsis a name.
    error : 1D-array
        If None, no errors will be plotted.
        If 1D array, error bars will be used for plotting.

    Yields
    ------

    A plot for each OGS folder in path_to_results.
    """

    import matplotlib.pyplot as plt
    import os.path

    if error == None:
        plt.plot(location, parameter, label=y_label)
    elif error != None:
        plt.errorbar(location, parameter, error, label=y_label)
    plt.xlabel("location [m]")
    plt.ylabel(y_label)
    plt.title(os.path.basename(path_to_results))
    plt.savefig(
        path_to_results + "/" + y_label + "_vs_location.png",
        bbox_inches="tight",
        dpi=300,
    )


def plot_error_vs_tc_t_and_s(
    results, path_to_results, lims=None, locs=None, comment="", abs=True, yaxis="log"
):
    """
    Plot errors of input and output parameters (S, T, tc) vs input tc.
    This results in 3 (S,T,tc) plots with multiple graphs according to the
    different observation points.

    Parameters
    ----------

    results : pandas dataframe
        Dataframe with results generated by multi_psd.py.
    path_to_results : string
        Path where to store the images
    abs : bool
        Take the mean over absolute errors or not.
    locs : list of integers
        default: [200, 400, 600, 800, 990]
    lims : list of two tuples with limits for x and y. [(x_min, x_max), (y_min, y_max)]
        Default: No limits
    yaxis : string, Default: "log"
        "log": loglog plot
        "lin": semilogx
    """

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import string

    from cycler import cycler
    # define a list for numbers for plots
    alpha_list = ["({})".format(i) for i in string.ascii_lowercase] + ["({})".format(i) for i in string.ascii_uppercase]
    
    try:
        if locs == None:
            locs = [200, 400, 600, 800, 990]
    except ValueError:
        pass

    # set matplotlib style
    CWD = os.getcwd()
    CWD_p = os.path.dirname(CWD)
    plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])

    fig, axes = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(15,7))

    if not os.path.exists(path_to_results + "/error_vs_tc"):
        os.mkdir(path_to_results + "/error_vs_tc")
    for i, (error, ax) in enumerate(zip(["err_S", "err_T"], axes.flat)):
        tc_agg = results["tc_in"].apply(lambda x: np.around(x, 2)).unique()
        tc_agg.sort()
        for loc in locs:  # results.obs_loc.unique():
            results_loc = results[error][results["obs_loc"] == loc]
            err_vs_tc_at_loc = []
            for tc in tc_agg:
                print(
                    "...currently grouping values for tc = " + str(tc),
                    "location: " + str(loc),
                    "Error: " + error,
                )
                # append error to list for specific tc
                if abs == True:
                    err_vs_tc_at_loc.append(
                        np.mean(
                            np.abs(
                                results_loc[
                                    results["tc_in"].apply(lambda x: np.around(x, 2))
                                    == tc
                                ]
                            )
                        )
                    )
                if abs == False:
                    err_vs_tc_at_loc.append(
                        np.mean(
                            results_loc[
                                results["tc_in"].apply(lambda x: np.around(x, 2)) == tc
                            ]
                        )
                    )
            if yaxis == "log":
                ax.loglog(tc_agg, err_vs_tc_at_loc, label=str(loc) + " m")
            elif yaxis == "lin":
                ax.semilogx(tc_agg, err_vs_tc_at_loc, label=str(loc) + " m")
            else:
                print("Please set argument 'yaxis' either to 'lin' or 'log'.")
        if lims != None:
            ax.set_ylim(lims[1])
            ax.set_xlim(lims[0])
        if error == "err_S":
            ax.set_ylabel(r"Error S [%]", fontsize=24)
        elif error == "err_T":
            ax.set_ylabel(r"Error T [%]", fontsize=24)
        ax.annotate(alpha_list[i], (10, 90), fontsize=16)
        ax.set_xlabel(r"$t_c\ [days]$", fontsize=24)
        ax.legend()
    plt.savefig(
        path_to_results
        + "/error_vs_tc/"
        + comment
        + "_"
        + error
        + "_vs_tc_abs_"
        + str(abs)
        + ".pdf",
        dpi=dpi,
        bbox_inches="tight",
    )
    print("Rdy.")
    plt.close()




if __name__ == "__main__":
    import pandas as pd
    import numpy as np


    # # Execute plot_error_vs_tc_t_and_s for spectral analysis paper 1
    # # --------------------------------------------------------------------------
    #results = pd.read_csv("/Users/houben/phd/studies/spectral_analysis/20190322_spectral_analysis_homogeneous_neu/combined_results/1_results_merge.csv")
    # path_to_results = "/Users/houben/phd/studies/spectral_analysis/20190322_spectral_analysis_homogeneous_neu/combined_results/error_vs_tc"
    # locs = [100, 200, 500, 800, 950]
    # #locs = [100]
    # lims = [(0.1, 10000), (0, 100)]
    # plot_error_vs_tc_t_and_s(results, path_to_results, lims=lims, locs=locs, comment="", abs=True, yaxis="lin")

    # check if the dataset is correct and the plot was correct 
    # # --------------------------------------------------------------------------
    # # --------------------------------------------------------------------------
    # results = pd.read_csv("/Users/houben/phd/studies/spectral_analysis/20190322_spectral_analysis_homogeneous_neu/combined_results/1_results_merge.csv")
    # path_to_results = "/Users/houben/phd/studies/spectral_analysis/20190322_spectral_analysis_homogeneous_neu/combined_results/error_vs_tc/NEW_TEST"
    # results.columns
    # # (['Unnamed: 0', 'Unnamed: 0.1', 'name', 'S_in', 'T_in', 'D_in', 'T_out',
    # #        'S_out', 'tc_out', 'cov', 'obs_loc', 'time_step_size', 'time_steps',
    # #        'model_period', 'which', 'recharge', 'aquifer_length',
    # #        'aquifer_thickness', 'D', 'D_cov', 'tc_in', 'err_S', 'err_T', 'err_tc'],
    # #       dtype='object')
    # from analysis.processing import percent_difference_fraction
    # a = percent_difference_fraction(results.T_in, results.T_out)
    # np.count_nonzero(np.around(results.err_T, 0) == np.around(a, 0))
    # results.err_T = percent_difference_fraction_log(results.T_in, results.T_out)
    # results.err_S = percent_difference_fraction_log(results.S_in, results.S_out)
    
    # from analysis.calc_tc import calc_tc   
    # # 'linear', 'gelhar', 'liang2013', 'dupuit' or 'liang2015'
    # tc_out_cal = calc_tc(1000, results.S_out, results.T_out, which="dupuit")
    # tc_in_cal = calc_tc(1000, results.S_in, results.T_in, which="dupuit")

    # # don't know how the tc_in was calculated, tc_out seems to be the correct tc
    # results.tc_in
    # results.tc_out
    
    # # found out that the name of the folders are wrong because they were
    # # calculated with the "gelhar" equation
    # [print(i) for i in results.name if i[:4] == "1313"]
    # tc_test = calc_tc(1000, 5.01e-05, 3.98e-04, which="gelhar")
    # locs = [100, 200, 500, 800, 950]
    # lims = [(0.1, 10000), (0, 100)]
    # plot_error_vs_tc_t_and_s(results, path_to_results, lims=lims, locs=locs, comment="", abs=True, yaxis="lin")

    # import matplotlib.pyplot as plt
    # plt.plot(tc_out_cal)
    # diff =  tc_out_cal - tc_in_cal
    # plt.plot(diff)
    # plt.plot(tc_in_cal)
    # plt.ylim(0,400000)
    # plt.plot(tc_in_cal)
    # plt.plot(results.tc_in)
    # plt.plot(tc_in_cal - results.tc_in)

    # plt.semilogy(sorted(results.tc_in))
    # plt.semilogy(sorted(tc_in_cal))
    # plt.show()
    # # --------------------------------------------------------------------------
    # # --------------------------------------------------------------------------



    # # Execute plot_error_vs_tc_t_and_s for spectral analysis paper 1 for the new mesh
    # # --------------------------------------------------------------------------
    results = pd.read_csv("/Users/houben/phd/studies/spectral_analysis/20210823_spectral_analysis_homogeneous_neu_new_mesh/combined_results_min/min_results_merge.csv")
    path_to_results = "/Users/houben/phd/studies/spectral_analysis/20210823_spectral_analysis_homogeneous_neu_new_mesh/combined_results_min"
    from analysis.calc_tc import calc_tc
    # 'linear', 'gelhar', 'liang2013', 'dupuit' or 'liang2015'
    # add a column for ic_in
    from analysis.processing import percent_difference_fraction, percent_difference_fraction_log
    results["tc_in"] = calc_tc(1000, results.S_in, results.T_in, which="liang2013")
    results["err_S"] = percent_difference_fraction(results.S_in, results.S_out)
    results["err_T"] = percent_difference_fraction(results.T_in, results.T_out)

    locs = [100, 200, 500, 800, 950]
    lims = [(0.01, 100000), (0, 100)]
    plot_error_vs_tc_t_and_s(results, path_to_results, lims=lims, locs=locs, comment="", abs=True, yaxis="lin")
        
    """
    # execute the plot_eff_thickness_vs_connection_length_and_thickness
    results = pd.read_csv("/Users/houben/phd/results/20191115_ogs_homogeneous_thickness/1_results_merge.csv")
    path = "/Users/houben/phd/results/20191115_ogs_homogeneous_thickness"
    obs_locs = [1000]
    stors = 0.01
    plot_eff_thickness_vs_connection_length_and_thickness(results, path, obs_locs, stors)
    """

    """
    # execute plot_kf_vs_connection_length_and_thickness(results, path)
    stors = [0.01, 0.0001]
    obs_locs = [0, 50, 100, 200, 500, 800, 900, 950, 1000]
    # obs_locs = [500]
    results = pd.read_csv(
        "/Users/houben/phd/results/20191115_ogs_homogeneous_thickness/1_results_merge.csv"
    )
    path = "/Users/houben/phd/results/20191115_ogs_homogeneous_thickness"
    plot_k_vs_connection_length_and_thickness(results, path, obs_locs, stors)
    """

    # 20190322_spectral_analysis_homogeneous_neu
    # --------------------------------------------
    # results = pd.read_csv("/Users/houben/phd/results/20190322_spectral_analysis_homogeneous_neu/combined_results/1_results_merge.csv")
    # path_to_results = "/Users/houben/phd/results/20190322_spectral_analysis_homogeneous_neu/combined_results"
    # locs = [100, 200, 500, 800, 950]
    # lims = [(0.1, 10000), (0, 100)]
    # plot_error_vs_tc(results, path_to_results, locs=locs, abs=True, lims=lims, yaxis="lin")
    #plot_errors_vs_loc(results, path_to_results)
    #plot_heatmap(results, path_to_results)

    # call plot_errors_vs_loc_aggregate for different parameter combinations
    # --------------------------------------------
    # bins for category "tc_in"
    # error = ["err_S", "err_T", "err_tc"]
    # aggregate = ["S_in", "T_in", "tc_in"]
    # bins = [
    # np.power(10,np.linspace(-5,-1,10)).tolist(),
    # np.power(10,np.linspace(-6,-2,10)).tolist(),
    # [1,10,100,1000,3000,9000,15000,20000,25000,30000,100000],
    # ]
    # for err in error:
    #    for i,agg in enumerate(aggregate):
    #        #print(err,agg,bins[i])
    #        plot_errors_vs_loc_aggregate(results, path_to_results, err, agg, bins[i], abs=True)


    ## execute plot_parameter_vs_location_block
    ## ----------------------------------------

    #borders = [800,810,820,830,840,850,860,870,880,890,900]
    #borders = [700,710,720,730,740,750,760,770,780,790,800]
    #borders = [200,210,250,280,290,300,310,320,330,340,500]
    #borders = [110,120,130,140,150,160,170,180,190,200]
    #borders = [50,100,200,300,400,500,600,800,950,970,990]
    #borders = [int(i) for i in np.linspace(10,990,99)]
    #borders = [10,20,30,40,50, 200, 400, 500,600, 700, 990]
    # borders = [50, 500, 990]
    #facecolor = "None"
    #recharge_list = ["recharge_daily.txt", "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # # # configuration for 20190531_SA_hetero_block
    # results = pd.read_csv("/Users/houben/phd/results/20190531_SA_hetero_block/results_merge.csv")
    # T_out_mhm = results.T_out[results.recharge == "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # path_to_results = "/Users/houben/phd/results/20190531_SA_hetero_block"
    # %matplotlib qt
    # plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, threshold=100, comment="S_0.003_border_", facecolor=facecolor)

    # configuration for 20190717_SA_hetero_block_2 for paper
    # --------------------------------------------
    #borders = [50,400,500,600,800,950,990]
    # borders = [50,500,990]
    # facecolor = "None"
    # recharge_list = ["recharge_daily.txt", "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # results = pd.read_csv("/Users/houben/phd/results/20190717_SA_hetero_block_2/1_results_merge.csv")
    # path_to_results = "/Users/houben/phd/results/20190717_SA_hetero_block_2"
    # plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="Paper_S_0.003_", facecolor=facecolor)

    # configuration for 20190531_SA_hetero_block for paper
    # --------------------------------------------
    #borders = [50,400,500,600,800,950,990]
    # borders = [50,500,990]
    # facecolor = "None"
    # recharge_list = ["recharge_daily.txt", "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # results = pd.read_csv("/Users/houben/phd/results/20190531_SA_hetero_block/results_merge.csv")
    # path_to_results = "/Users/houben/phd/results/20190531_SA_hetero_block"
    # plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="Paper_S_0.003_", facecolor=facecolor)

    # configuration for 20190531_SA_hetero_block for Supplement of paper
    # --------------------------------------------
    # borders = [50,100,200,300,400,500,600,800,950,970,990]
    # facecolor = "None"
    # recharge_list = ["recharge_daily.txt", "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # results = pd.read_csv("/Users/houben/phd/studies/spectral_analysis/20190531_SA_hetero_block/results_merge.csv")
    # path_to_results = "/Users/houben/phd/studies/spectral_analysis/20190531_SA_hetero_block"
    # plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="Paper_Supp_S_0.003_", facecolor=facecolor)

    # configuration for 20190717_SA_hetero_block_2 for paper
    # --------------------------------------------
    # borders = [50,100,200,300,400,500,600,800,950,970,990]
    # facecolor = "None"
    # recharge_list = ["recharge_daily.txt", "recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    # results = pd.read_csv("/Users/houben/phd/studies/spectral_analysis/20190717_SA_hetero_block_2/1_results_merge.csv")
    # path_to_results = "/Users/houben/phd/studies/spectral_analysis/20190717_SA_hetero_block_2"
    # plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="Paper_Supp_S_0.003_", facecolor=facecolor)

    # for i in range(0,10):
    #     intervall = i * 100
    #     start = 10 + intervall
    #     end = 90 + intervall
    #     borders = [int(i) for i in np.linspace(start,end,9)]
    #     print(borders)
    #     plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="S_0.003_border_"+str(start)+"_"+str(end)+"_", facecolor=facecolor)
    #     #plot_parameter_vs_location_block(results, path_to_results, borders=borders, S_in = 0.003, recharge_list=recharge_list, comment="S_0.003_border_"+str(start)+"_"+str(end)+"_", facecolor=facecolor)

    # execute plot_baseflow_sa_vs_boundary_block
    # --------------------------------------------
    # import pandas as pd
    # results = pd.read_csv("/Users/houben/phd/results/20190717_SA_hetero_block_2/combined_results/baseflow_results_merge.csv")
    # plot_baseflow_sa_vs_boundary_block(results,"/Users/houben/phd/results/20190717_SA_hetero_block_2/combined_results/baseflow_results_merge.csv")
