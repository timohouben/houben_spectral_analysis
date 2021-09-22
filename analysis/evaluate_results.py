import numpy as np
dpi = 250

def identify_misfits(path, filename="results.csv"):
    """
    Identifies the model runs which have higher covariance matrix than 1 and prints input parameters.

    Parameters
    ----------

    path : string
        Path to results
    filename : string
        "results.csv" is default
    """
    count = 0
    import pandas as pd

    df = pd.read_csv(path + "/" + filename)
    for i, item in enumerate(df["cov"][:]):
        try:
            if float(item[3:17]) > 1:
                print(
                    item[3:17],
                    i,
                    df["name"][i],
                    "S_in: ",
                    df["S_in"][i],
                    "T_in: ",
                    df["T_in"][i],
                )
                count = count + 1
        except TypeError:
            pass
        except ValueError:
            pass
    print(str(count), "misfits identified")


def evaluate_fit(path_to_results, filename="results.csv"):
    """
    Evaluate the fit of power spectra.

    Parameters
    ----------
    path : string
        Path to results and where to store the plots.
    filename : string
        Name of file. Default "results.csv".
    """
    import pandas as pd

    results = pd.read_csv(path_to_results + "/" + filename)

    def plot(pivotted, variance):
        import seaborn as sns
        import os
        import numpy as np
        from matplotlib.colors import LogNorm
        import math

        barmin, barmax = 1e-18, 1e-8
        cbar_ticks = [1e-20, 1e-18, 1e-16, 1e-14, 1e-12, 1e-10]
        log_norm = LogNorm(vmin=barmin, vmax=barmax)
        ax = sns.heatmap(
            pivotted,
            cmap="coolwarm",
            vmax=barmax,
            vmin=barmin,
            norm=log_norm,
            cbar_kws={"ticks": cbar_ticks},
        )  # , yticklabels=achsislabel_y, xticklabels=achsislabel_x)
        # ax.invert_yaxis()
        fig = ax.get_figure()
        if not os.path.exists(path_to_results + "/heatmap_variance"):
            os.mkdir(path_to_results + "/heatmap_variance")

        fig.savefig(
            path_to_results + "/heatmap_variance" + "/" + str(obs_loc) + "_" + variance,
            dpi=dpi,
        )
        fig.clf()

    from processing import identify_numbers_from_string

    for obs_loc in results["obs_loc"].unique():
        # extract only rows with obs_loc==obs_loc
        df_obs_loc = results[results.obs_loc == obs_loc]
        # extract columns for plotting
        df_obs_loc_cut = df_obs_loc[["S_in", "T_in", "cov"]]
        # get values for sigma S and sigma T seperately from column cov
        df_obs_loc_cut["cov_numbers"] = df_obs_loc_cut["cov"].apply(
            identify_numbers_from_string
        )
        df_obs_loc_cut["sigma_S"] = df_obs_loc_cut["cov_numbers"].apply(lambda x: x[0])
        df_obs_loc_cut["sigma_T"] = df_obs_loc_cut["cov_numbers"].apply(lambda x: x[3])
        # convert objects to floats
        df_obs_loc_cut.sigma_S = pd.to_numeric(df_obs_loc_cut.sigma_S)
        df_obs_loc_cut.sigma_T = pd.to_numeric(df_obs_loc_cut.sigma_T)
        for variance in ["sigma_S", "sigma_T"]:
            pivot_df_obs_loc_cut = df_obs_loc_cut.pivot("S_in", "T_in", variance)
            # plot heatmap
            import numpy as np

            plot(pivot_df_obs_loc_cut, variance)


if __name__ == "__main__":
    pass
    #evaluate_fit("/Users/houben/Desktop/test", "csv_merge.csv")
