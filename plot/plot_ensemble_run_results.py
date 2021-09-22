def plot_combined_violins_for_ensemble_runs_hetero(
    path, filename, order, len_scales, y, savename="", y_lims=None, bw=None, yachsis="lin", format=".png"
):
    """
    Make violin plots for ensemble runs with log-normal distributed hydraulic
    conductivity. Each side of a violin represents a
    correlation length. Each violin represents an observation point.

    Parameters
    ----------

    path : string
        Path to results.
    filename : string
        Name of the results file.
    savename : string
        Specify a name for the file to be saved
    order : list
        Give a list with the observation points
    len_scales : list
        give a list of size 2 with correlation lenghts
    y : string
        Which column of data frame to use for violins.
    ylims : tuple (y_min,y_max), Default: None
        Tuple of limits for y-axis:
    bw : float, Default: None
        Bandwith of violin plots.
    yachsi : string, Default: "lin"
        "lin" : linear y-achsis
        "log" : log y-achsis


    Yields
    ------
    A violin plot in the provided directory.

    """

    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    import numpy as np

    # ------------------------
    # matplotlib costum design
    # ------------------------
     # set matplotlib style
    CWD = os.getcwd()
    CWD_p = os.path.dirname(CWD)
    plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
    # ------------------------
    # plt.style.use("ggplot")
    # plt.style.available
    # # plt.rc('text', usetex=True)
    # plt.rc("font", family="serif")
    # plt.rc("xtick", labelsize="x-large")
    # plt.rc("ytick", labelsize="x-large")
    # fig = plt.figure(figsize=(16, 9))

    data = pd.read_csv(path + "/" + filename)
    fig = plt.figure()
    sns.set_context("paper")
    ax1 = sns.catplot(
        x="obs_loc",
        y=y,
        order=order,
        data=data[data.len_scale.isin(len_scales)],
        kind="violin",
        height=5,
        aspect=2,
        scale="count",
        hue="len_scale",
        split=True,
        legend=False,
        bw=bw,
    )
    if y == "T_out":
        ax2 = plt.hlines(
            y=data.T_in_geo.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="geo. mean",
            color="c",
            alpha=0.5,
        )
        ax2 = plt.hlines(
            y=data.T_in_har.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="har. mean",
            color="r",
            alpha=0.005,
        )
        ax2 = plt.hlines(
            y=data.T_in_ari.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="ari. mean",
            color="y",
            alpha=0.005,
        )
        ax1.set_ylabels(r"T $\rm[m^2s^{-1}]$", fontsize=24)
    if y == "S_out":
        ax2 = plt.hlines(
            y=data.S_in.unique(),
            xmin=-0.5,
            xmax=4.5,
            label="S input",
            color="c",
            alpha=0.5,
        )
        ax1.set_ylabels(r"S $\rm[-]$", fontsize=24)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    leg = plt.legend(loc="lower left", fontsize="large")
    for line in leg.get_lines():
        line.set_alpha(1)
    ax1.set_xlabels("location of observation point [m]", fontsize=24)
    if yachsis == "log":
        plt.yscale('log')
    elif yachsis == "lin":
        pass
    else:
        print("Argument yachsis has to be 'lin' or 'log'. Selected default 'lin'.")
    plt.ylim(y_lims)
    plt.savefig(
        path
        + "/"
        + savename
        + y
        + "_len_"
        + "_".join([str(i) for i in len_scales])
        + "_obs_"
        + "_".join([str(i) for i in order])
        + "_bw_"
        + str(bw)
        + "_"
        + yachsis
        + format,
        dpi=300,
        bbox_inches="tight",
    )
    print("The ari/geo/har means are:")
    print(np.mean(data.T_in_ari.unique()))
    print(np.mean(data.T_in_geo.unique()))
    print(np.mean(data.T_in_har.unique()))



def plot_combined_violins_for_ensemble_runs_layered(
    path, filename, order, recharges, y, S_in, savename="", y_lims=None, ylog=True, bw=None
):
    """
    Make violin plots for ensemble runs with layered aquifers.
    Each side of a violin represents a recharge process.
    Each violin represents an observation point.

    Parameters
    ----------

    path : string
        Path to results.
    filename : string
        Name of the results file.
    savename : string
        Specify a name for the file to be saved
    order : list
        Give a list with the observation points
    recharges : list
        give a list of size 2 with names of recharge-files
    y : string
        Which column of data frame to use for violins.
    S_in : float
        For which storage value shall data be plotted?
    ylims : tuple (y_min,y_max)
        Tuple of limits for y-axis:
    ylog : Bool (default: True)
        If True, log scale for y-achsis is used
    bw : float
        Bandwith of violin plots.


    Yields
    ------
    A violin plot in the provided directory.

    """

    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt

    data = pd.read_csv(path + "/" + filename)
    data = data[data["S_in"] == S_in]
    fig = plt.figure()
    sns.set_context("paper")
    ax1 = sns.catplot(
        x="obs_loc",
        y=y,
        order=order,
        data=data[data.recharge.isin(recharges)],
        kind="violin",
        height=5,
        aspect=2,
        scale="count",
        hue="recharge",
        split=True,
        legend=False,
        bw=bw,
    )
    if y == "T_out":
        ax2 = plt.hlines(
            y=data.T_in_geo.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="geomean",
            color="c",
            alpha=0.5,
        )
        ax3 = plt.hlines(
            y=data.T_in_har.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="harmean",
            color="r",
            alpha=0.5,
        )
        ax4 = plt.hlines(
            y=data.T_in_ari.unique(),
            xmin=-0.5,
            xmax=len(order) - 0.5,
            label="arimean",
            color="y",
            alpha=0.5,
        )
        ax1.set_ylabels("$T\;[m^2/s]$ derived by fit")
    if y == "S_out":
        ax2 = plt.hlines(
            y=data.S_in.unique(),
            xmin=-0.5,
            xmax=4.5,
            label="S input",
            color="c",
            alpha=0.5,
        )
        ax1.set_ylabels("$S\;[-]$ derived by fit")
    plt.legend(loc="upper left")
    ax1.set_xlabels("location of observation point [m]")
    plt.ylim(y_lims)
    plt.title("Derived Transmissivity for Storativity " + str(S_in))
    if ylog == True:
        plt.yscale("log")
    plt.savefig(
        path
        + "/"
        + savename
        + y
        + "_S_"
        + str(S_in)
        + "_len_"
        + "_".join([str(i) for i in recharges])
        + "_obs_"
        + "_".join([str(i) for i in order])
        + "_ylog_"
        + str(ylog)
        + "_bw_"
        + str(bw)
        + ".png",
        dpi=300,
        bbox_inches="tight",
    )


def plot_violins_for_ensemble_runs_layered_baseflow(path, filename, recharges, y, savename="", y_lims=None, ylog=True, bw=None):
    """
    Plot violin plot for each setup (e.g. whitenoise and stor == 0.3)
    next to each other. Categorial violin plot.

    Parameters
    ----------

    path : string
        Path to results.
    filename : string
        Name of the results file.
    savename : string
        Specify a name for the file to be saved
    order : list
        Give a list with the observation points
    recharges : list
        give a list of size 2 with names of recharge-files
    y : string
        Which column of data frame to use for violins.
    S_in : float
        For which storage value shall data be plotted?
    ylims : tuple (y_min,y_max)
        Tuple of limits for y-axis:
    ylog : Bool (default: True)
        If True, log scale for y-achsis is used
    bw : float
        Bandwith of violin plots.

    Yields
    ------
    A violin plot in the provided directory.

    """

    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt

    data = pd.read_csv(path + "/" + filename)
    fig = plt.figure()
    sns.set_context("paper")
    ax1 = sns.catplot(
        x="S_in",
        y=y,
        order=data.S_in.unique(),
        data=data[data.recharge.isin(recharges)],
        kind="violin",
        height=5,
        aspect=2,
        scale="count",
        hue="recharge",
        split=True,
        legend=False,
        bw=bw,
    )
    if y == "T_out":
        ax2 = plt.hlines(
            y=data.T_in_geo.unique(),
            xmin=-0.5,
            xmax=2 - 0.5,
            label="geomean",
            color="c",
            alpha=0.5,
        )
        ax3 = plt.hlines(
            y=data.T_in_har.unique(),
            xmin=-0.5,
            xmax=2 - 0.5,
            label="harmean",
            color="r",
            alpha=0.5,
        )
        ax4 = plt.hlines(
            y=data.T_in_ari.unique(),
            xmin=-0.5,
            xmax=2 - 0.5,
            label="arimean",
            color="y",
            alpha=0.5,
        )
        ax1.set_ylabels("$T\;[m^2/s]$ derived by fit")
    if y == "S_out":
        ax2 = plt.hlines(
            y=data.S_in.unique(),
            xmin=-0.5,
            xmax=4.5,
            label="S input",
            color="c",
            alpha=0.5,
        )
        ax1.set_ylabels("$S\;[-]$ derived by fit")
    plt.legend(loc="upper left")
    ax1.set_xlabels("input Storativity [-]")
    plt.ylim(y_lims)
    plt.title("Derived Transmissivity from Baseflow")
    if ylog == True:
        plt.yscale("log")
    plt.savefig(
        path
        + "/"
        + savename
        + y
        + "_baseflow_"
        + "_len_"
        + "_".join([str(i) for i in recharges])
        + "_ylog_"
        + str(ylog)
        + "_bw_"
        + str(bw)
        + ".png",
        dpi=300,
        bbox_inches="tight",
    )

def plot_combined_violins_for_ensemble_runs_hetero_subplots(
    path, filename, order, len_scales, y, savename="", y_lims=None, bw=None, yachsis="lin", format=".pdf"
):
    """
    Make violin plots for ensemble runs with log-normal distributed hydraulic
    conductivity. Each side of a violin represents a
    correlation length. Each violin represents an observation point.

    THIS FUNCTIONCREATES TO MANY PLOTS FOR ANY REASON!!
    I DIDN'T FIGURE OUT WHY?!?!? Saved PDF from interactive view and edited with 
    incscape.

    Parameters
    ----------

    path : string
        Path to results.
    filename : string
        Name of the results file.
    savename : string
        Specify a name for the file to be saved
    order : list
        Give a list with the observation points
    len_scales : list
        give a list of size 2 with correlation lenghts
    ys : list of strings
        Which column of data frame to use for violins.
    ylims : tuple (y_min,y_max), Default: None
        Tuple of limits for y-axis:
    bw : float, Default: None
        Bandwith of violin plots.
    yachsi : string, Default: "lin"
        "lin" : linear y-achsis
        "log" : log y-achsis


    Yields
    ------
    A violin plot in the provided directory.

    """

    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    import string

    # ------------------------
    # set matplotlib style
    CWD = os.getcwd()
    CWD_p = os.path.dirname(CWD)
    plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
    # define a list for numbers for plots
    alpha_list = ["({})".format(i) for i in string.ascii_lowercase] + ["({})".format(i) for i in string.ascii_uppercase]
    # fig = plt.figure(figsize=(16, 9))
    # ------------------------

    data = pd.read_csv(path + "/" + filename)
    #fig = plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(15, 15))
    for i, (ax, y) in enumerate(zip(axes, ys)):

        ax_temp = sns.catplot(
            x="obs_loc",
            y=y,
            order=order,
            data=data[data.len_scale.isin(len_scales)],
            kind="violin",
            height=5,
            aspect=2,
            scale="count",
            hue="len_scale",
            split=True,
            legend=False,
            bw=bw,
            ax=ax,
        )
        # ax1 = plt.subplot(211)
        # ax1 = ax_temp
        # print(ax_temp)
        print(fig)
        #plt.delaxes(a)
        if y == "T_out":
            ax.hlines(
                y=data.T_in_geo.unique(),
                xmin=-0.5,
                xmax=len(order) - 0.5,
                label="geo. mean",
                color="c",
                alpha=0.5,
            )
            ax.hlines(
                y=data.T_in_har.unique(),
                xmin=-0.5,
                xmax=len(order) - 0.5,
                label="har. mean",
                color="r",
                alpha=0.005,
            )
            ax.hlines(
                y=data.T_in_ari.unique(),
                xmin=-0.5,
                xmax=len(order) - 0.5,
                label="ari. mean",
                color="y",
                alpha=0.005,
            )
            ax.set_ylabel(r"$\rm T\;[m^2s^{-1}]$", fontsize=24)
            ax.annotate(alpha_list[i], (0, 7e-3), fontsize=16)
        if y == "S_out":
            ax.hlines(
                y=data.S_in.unique(),
                xmin=-0.5,
                xmax=4.5,
                label="S input",
                color="c",
                alpha=0.5,
            )
            ax.set_ylabel(r"$\rm S\;[-]$", fontsize=24)
            ax.set_xlabel("")
            # #ax.get_xaxis().set_visible(False)
            ax.annotate(alpha_list[i], (0, 1.8e-1), fontsize=16)
        if yachsis == "log":
            ax.set_yscale('log')
        elif yachsis == "lin":
            pass
        else:
            print("Argument yachsis has to be 'lin' or 'log'. Selected default 'lin'.")
        leg = ax.legend(loc="lower left")
        for line in leg.get_lines():
            line.set_alpha(1)
    plt.subplots_adjust(hspace=0.01)
    ax.set_xlabel("location of observation point [m]", fontsize=24)
    #ax.set_ylim(y_lims)
    #plt.show()
    plt.savefig(
        path
        + "/"
        + savename
        + "_".join(ys)
        + "_len_"
        + "_".join([str(i) for i in len_scales])
        + "_obs_"
        + "_".join([str(i) for i in order])
        + "_bw_"
        + str(bw)
        + "_"
        + yachsis
        + format,
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == '__main__':
    pass
    #    order = [[],[]]
    #    len_scale = [[],[]]
    #    y = [,]
    #    for orderT in order:
    #        for len_scalesT in len_scales:
    #            for yT in y:
    #                plot_combined_violins_for_ensemble_runs("/Users/houben/phd/results/20190513_spec_anal_hetero_ensemble_1", "merge_results.csv", orderT, len_scalesT,yT,bw=0.2)


    # plot the results for 20191108_results_merge
    # and for spectral analysis paper 1
    # coded the plot_combined_violins_for_ensemble_runs_hetero_subplots
    # to plot both in one plot but didn't work so went back to
    # plot_combined_violins_for_ensemble_runs_hetero
    # and plotted S and T seperately
    # ------------------------------------------------------------------------------
    # path = "/Users/houben/phd/studies/spectral_analysis/20190513_spec_anal_hetero_ensemble_1/20191114_combined_results"
    # filename = "20191108_results_merge.csv"
    # order = [100,200,500,800,940]
    # len_scales = [5, 15]
    # ys = ["S_out", "T_out"]
    # yachsis = "log"

    # plot_combined_violins_for_ensemble_runs_hetero_subplots(path, filename, order, len_scales, ys, savename="TEST", y_lims=None, bw=None, yachsis=yachsis, format=".pdf")


    # plot the results for 20191108_results_merge
    path = "/Users/houben/phd/studies/spectral_analysis/20190513_spec_anal_hetero_ensemble_1/20191114_combined_results"
    filename = "20191108_results_merge.csv"
    order = [100,200,500,800,940]
    len_scales = [5, 15]
    y = "T_out"
    yachsis = "log"


    plot_combined_violins_for_ensemble_runs_hetero(path, filename, order, len_scales, y, savename="TEST", y_lims=None, bw=None, yachsis=yachsis, format=".pdf")




    '''
    # plot the results for 20190917_ogs_layered_ensemble

    path = "/Users/houben/phd/results/20190917_ogs_layered_ensemble"
    filename = "1_results_merge.csv"
    order = [0,50,250,500,750,900,990]
    recharges = ["recharge_daily.txt","recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    y = "T_out"
    savename="plot"
    # y lims for log plot incl all means
    #y_lims=[7e-5,1e-1]
    # y lims for log plot zoomed in
    #y_lims=[3e-3,8e-3]
    # y lims for lin plot
    y_lims=[0.00002,0.008]
    S_in = 0.003
    ylog=False

    plot_combined_violins_for_ensemble_runs_layered(path, filename, order, recharges, y, S_in=S_in, savename="", y_lims=y_lims, ylog=ylog, bw=None)
    '''
    '''
    # plt the baseflow results for 20190917_ogs_layered_ensemble

    path = "/Users/houben/phd/results/20190917_ogs_layered_ensemble"
    filename = "1_results_merge.csv"
    recharges = ["recharge_daily.txt","recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt"]
    y = "T_out"
    # lims for lin axis
    #y_lims = [0.00002,0.008]
    # lims for log axis
    #y_lims=[7e-5,1e-1]
    # lims for log axis zoomed
    y_lims=[3e-3,8e-3]
    plot_violins_for_ensemble_runs_layered_baseflow(path, filename, recharges, y, savename="", y_lims=y_lims, ylog=True, bw=None)
    '''
