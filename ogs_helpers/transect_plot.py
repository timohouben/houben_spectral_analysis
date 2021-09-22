def extract_timeseries(path, which="mean", process="GROUNDWATER_FLOW"):
    """
    Function to extract time series for each observation point and store them as .txt.

    Parameters
    ----------

    path : string
        path to ogs directory
    which : string, 'mean', 'max', 'min'
        Which value should be taken from the vertical polyline of obs point.
    process : string
        Which modelling process.

    Yields
    ------

    .txt from every observation point with head values

    """

    import numpy as np
    from ogs5py.reader import readtec_polyline
    import glob
    import vtk

    from ogs_helpers.tools import get_ogs_task_id

    # pipe vtk output errors to file
    errOut = vtk.vtkFileOutputWindow()
    errOut.SetFileName(path + "/VTK_Error_Out.txt")
    vtkStdErrOut = vtk.vtkOutputWindow()
    vtkStdErrOut.SetInstance(errOut)

    task_id = get_ogs_task_id(path)

    # read all tec files
    print("Reading tec-files from " + path)
    # This will throw an ERROR which is redirected to VTK_Error_Out.txt
    tecs = readtec_polyline(task_id=task_id, task_root=path)

    # extract the time series and save them as .txt
    for obs in tecs["GROUNDWATER_FLOW"].keys():
        time_steps = len(tecs["GROUNDWATER_FLOW"][obs]["TIME"])
        number_of_columns = tecs[process][obs]["HEAD"].shape[1]
        if which == "max":
            # select the maximum value (i.e. the uppermost) of polyline as long as polylines are defined from bottom to top
            head_ogs_timeseries_each_obs = tecs[process][obs]["HEAD"][
                :, number_of_columns - 1
            ]
        elif which == "min":
            # select the minimum value (i.e. the lowermost) of polyline as long as polylines are defined from bottom to top
            head_ogs_timeseries_each_obs = tecs[process][obs]["HEAD"][:, 0]
        elif which == "mean":
            head_ogs_timeseries_each_obs = []
            for step in range(time_steps):
                # calculates the mean of each time step
                head_ogs_timeseries_each_obs.append(
                    np.mean(tecs[process][obs]["HEAD"][step, :])
                )
            head_ogs_timeseries_each_obs = np.asarray(head_ogs_timeseries_each_obs)
        np.savetxt(
            str(path) + "/" + "head_ogs_" + str(obs) + "_" + str(which) + ".txt",
            head_ogs_timeseries_each_obs,
        )


def extract_rfd(path, rfd=1, export=True):
    """
    Function to extract the x and y values of a given rfd or to load it from previosly exportet .txt-files.

    Parameters
    ----------

    path : string
        path to ogs directory
    rfd : int
        Extract the #rfd curve from the rfd file and save it as txt. rfd = 0 : no extraction
    export : bool (Default: True)
        If True .txt-files will be saved to the path.

    Yields
    ------

    A txt file for x and y values of given rfd and returns x_values and y_values as tuple of lists.

    """

    import numpy as np
    import os.path
    from ogs_helpers.tools import get_ogs_task_id

    # get the task if from ogs model run
    task_id = get_ogs_task_id(path)

    # extract the series given in the rfd file
    if rfd != 0:
        if os.path.exists(str(path) + "/" + "rfd_curve#" + str(rfd) + "_x_values.txt") and os.path.exists(str(path) + "/" + "rfd_curve#" + str(rfd) + "_y_values.txt"):
            print("Txt.file of corresponding rfd curve already exists. Continuing without checking if content is correct.")
            # load this files to return from function
            x_values = np.loadtxt(path + "/" + "rfd_curve#" + str(rfd) + "_x_values.txt")
            y_values = np.loadtxt(path + "/" + "rfd_curve#" + str(rfd) + "_y_values.txt")
        else:
            from ogs5py import OGS
            ogs = OGS(task_root=path + "/", task_id=task_id)
            ogs.rfd.read_file(path=path + "/" + task_id + ".rfd")
            #print(ogs.rfd.get_block(rfd-1)[''])
            y_values = np.asarray([y_values[1] for y_values in ogs.rfd.get_block(rfd - 1)[""]])
            x_values = np.asarray([x_values[0] for x_values in ogs.rfd.get_block(rfd - 1)[""]])
            if export == True:
                np.savetxt(str(path) + "/" + "rfd_curve#" + str(rfd) + "_x_values.txt", x_values)
                np.savetxt(str(path) + "/" + "rfd_curve#" + str(rfd) + "_y_values.txt", y_values)
            else:
                pass
    return x_values, y_values



def plot_head_timeseries_vs_recharge(path, which="mean"):
    """
    Plot the head.txt with recharge time series on 2nd axis.
    Assuming that recharge is given in m/s and the time steps are in seconds with a dayli increment.
    """

    import glob
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import os.path

    # glob the name of the .txt files
    file_names_obs_list = glob.glob(path + "/*obs*" + which + "*.txt")
    if file_names_obs_list == []:
        print("No head.txt files have been extracted previously. Running function extract_timeseries with default parameters.")
        extract_timeseries(path, which="mean", process="GROUNDWATER_FLOW")
        file_names_obs_list = glob.glob(path + "/*obs*" + which + "*.txt")

    file_names_obs_list.sort()

    try:
        file_name_values = str(glob.glob(path + "/*1_y_values.txt")[0])
        file_name_time = str(glob.glob(path + "/*1_x_values.txt")[0])
    except IndexError:
        print("No rfd.txt files have been extracted previously. Running extract_rfd with default parameters.")
        extract_rfd(path, rfd=1, export=True)
        file_name_values = str(glob.glob(path + "/*1_y_values.txt")[0])
        file_name_time = str(glob.glob(path + "/*1_x_values.txt")[0])

    rfd = np.loadtxt(file_name_values)
    time = np.loadtxt(file_name_time)

    # calculate mm/day from m/s
    rfd = rfd * 86400 * 1000
    time = time / 86400

    # first axis for recharge
    fig, ax1 = plt.subplots(figsize=(14, 10))
    plt.title("head timeseries at different observations points")
    color = "tab:blue"
    ax1.set_xlabel("time [day]")
    ax1.set_ylabel(
        "recharge [mm/day]", color=color
    )  # we already handled the x-label with ax1
    ax1.bar(time, rfd, width=0.8, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    # ax1.set_ylim(0, 2.5)
    # ax1.set_ylim([min(recharge), max(recharge)*2])
    # ax1.set_yticks([0, 1, 2, 3, 4, 5])

    # second axis for head
    ax2 = ax1.twinx()
    # ax2.set_ylim(30,30.1)
    # ax2.set_yticks(np.arange(26,40,0.5))
    color = "tab:red"

    print("plotting...")
    if which == "min":
        for obs in file_names_obs_list:
            # derive the head for the given observation point from ogs
            head_ogs = np.loadtxt(obs)
            ax2.plot(time[:len(head_ogs)], head_ogs, label=str(obs)[-13:-5] + " OGS", linestyle="-")
    elif which == "max":
        for obs in file_names_obs_list:
            # derive the head for the given observation point from ogs
            head_ogs = np.loadtxt(obs)
            ax2.plot(time[:len(head_ogs)], head_ogs, label=str(obs)[-13:-5] + " OGS", linestyle="-")
    elif which == "mean":
        for obs in file_names_obs_list:
            # derive the head for the given observation point from ogs
            head_ogs = np.loadtxt(obs)
            ax2.plot(time[:len(head_ogs)], head_ogs, label=str(obs)[-14:-4] + " OGS", linestyle="-")

    ax2.set_ylabel("head [m]", color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    ax2.grid(color="grey", linestyle="--", linewidth=0.5, which="both")
    handles, labels = ax2.get_legend_handles_labels()
    ax1.legend(handles, labels, loc=6, facecolor="white", framealpha=100)

    fig.tight_layout()
    #plt.show()
    # make a string from list obs_per_plot
    fig.savefig(
        str(path)
        + "/"
        + str(os.path.basename(str(path)))
        + "_"
        + str(which)
        + "_"
        + ".png"
    )
    plt.close("all")
