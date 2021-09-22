from ogs5py.reader import readtec_polyline
import numpy as np
import os
import matplotlib.pyplot as plt
from transect_plot import extract_rfd
import sys

def get_start_time(rfd):
    """
    Identifies time of first increase in y_values of rfd.

    Parameters
    ----------

    rfd : 2D-array shape(x,2)
        [:,0] : times
        [:,1] : corresponding values
    """
    # generate a test
    # rfd = np.asarray([[i for i in range(0,100)],[10]*50+[20]*50]).transpose()
    # load and test
    # times_rfd, concentration_in = extract_rfd("/Users/houben/phd/modelling/20190602_Tc_vs_Tr/transport/setup/2b_transport_stor_0.001_kf_0.001_disp_10_1_rech_1.16e-08_time_1000_8640000", rfd=1, export=True)
    # rfd = np.asarray([times_rfd, concentration_in]).transpose()

    for i in range(0, len(rfd)):
        if i < len(rfd):
            if rfd[i, 1] < rfd[i + 1, 1]:
                input_time = rfd[i + 1, 0]
                break
        else:
            print("No start time identified.")
    return input_time


def get_avg_conc_and_times_from_polyline(task_id, task_root, single_file):
    """
    Parameters
    ----------
    task_id
    task_root
    single_file

    Yields
    ------
    avg_concentration : time series of averaged concentration for a polyline
    times : coresponding times
    time_steps : number of time steps
    nodes : number of nodes along the polyline

    """
    tec = readtec_polyline(
        task_id=task_id, task_root=task_root, single_file=single_file
    )

    time_steps = tec[var_name].shape[0]
    nodes = tec[var_name].shape[1]
    times = tec["TIME"]

    avg_concentration = []

    for i in range(0, time_steps):
        # get the node values of the velocities for ith time step
        node_concentration = tec[var_name][i, :]
        avg_concentration.append(np.mean(node_concentration))

    return avg_concentration, times, time_steps, nodes


def scenario_2d(task_id, task_root, single_file, var_name, normalize):
    avg_concentration, times, time_steps, nodes = get_avg_conc_and_times_from_polyline(
        task_id=task_id, task_root=task_root, single_file=single_file
    )
    times_rfd, concentration_in = extract_rfd(task_root, rfd=1, export=True)

    plt.figure(figsize=(10, 7))
    plt.plot(
        times_rfd / normalize,
        concentration_in,
        linewidth=4,
        label="input concentration",
        color="blue",
    )
    plt.plot(
        times / normalize,
        avg_concentration,
        linewidth=2,
        label="averaged outflow concentration",
        color="orange",
    )
    plt.xlabel("days")
    plt.ylabel("concentration")
    plt.title(os.path.basename(task_root))
    plt.legend()
    plt.savefig(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.png",
        dpi=300,
    )
    #plt.show()
    np.savetxt(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.txt",
        avg_concentration,
    )
    return avg_concentration


def scenario_2c(
    task_id,
    task_root,
    single_file,
    var_name,
    normalize,
    percentage,
    ref_conc=None,
):
    """
    Script to plot the averaged concentration along a polyline vs the input concentration.

    Parameters
    ----------

    ref_conc : float, default: None
        If None, its extrated from second block in .ic from ogs model run.

    """

    if ref_conc == None:
        from ogs5py import OGS
        import glob

        ogs = OGS(
            task_root=task_root,
            task_id=os.path.basename(str(glob.glob(task_root + "/*.bc")))[:-5],
            output_dir=task_root,
        )
        ogs.ic.read_file(path=task_root + "/" + task_id + ".ic")
        ref_conc = ogs.ic.get_block(1)["DIS_TYPE"][0][1]

    avg_concentration, times, time_steps, nodes = get_avg_conc_and_times_from_polyline(
        task_id=task_id, task_root=task_root, single_file=single_file
    )
    times_rfd, concentration_in = extract_rfd(task_root, rfd=1, export=True)
    start_time = get_start_time(np.asarray([times_rfd, concentration_in]).transpose())

    # 95 % of concentration
    for time, conc_out in zip(times, avg_concentration):
        if percentage * conc_out <= ref_conc and time > start_time:
            delta_t = time - start_time
            stop_time = delta_t + start_time
            print(
                "Avegerage concentration at the outlet has reached "
                + str(percentage)
                + " % of reference concentration "
                + str(ref_conc)
                + " after "
                + str(delta_t / normalize)
                + " days."
            )
            break
        if time == np.max(times):
            print(
                "Concentration has not reached "
                + str(percentage)
                + " % of reference concentration "
                + str(ref_conc)
                + " within modelling period. Reached {:.3f} %.".format(
                    ref_conc / conc_out
                )
            )
            delta_t = " greater than modelling period..."

    fig, ax1 = plt.subplots(figsize=(10, 7))
    ax1.set_title(os.path.basename(task_root))
    ax1.set_xlabel("days")
    ax1.set_ylabel("concentration in", color="b")
    ax1.tick_params("y", colors="b")
    ax1.plot(
        times / normalize,
        concentration_in,
        linewidth=4,
        label="input concentration",
        color="b",
    )
    ax1.legend(loc=2)
    if delta_t != " greater than modelling period...":
        plt.vlines(
            x=[start_time / normalize, stop_time / normalize],
            ymin=0,
            ymax=np.max(concentration_in),
            color="g",
        )
        plt.annotate(
            "  $\Delta t$ = "
            + str(delta_t / normalize)
            + " to reach "
            + str(percentage)
            + " % of reference concentration.",
            (stop_time / normalize, 0.3 * np.max(concentration_in)),
        )
    ax2 = ax1.twinx()
    ax2.plot(
        times / normalize,
        avg_concentration,
        linewidth=2,
        label="averaged outflow concentration",
        color="orange",
    )
    ax2.set_ylabel("concentration out", color="orange")
    ax2.tick_params("y", colors="orange")

    # if delta_t == " greater than modelling period...":
    #    plt.annotate(str(percentage) + " % is" + delta_t + " \nReached {:.3f} % in modelling period.".format(conc_out / conc_in), (np.median(times) + 1 ,0.8*max(concentration_in)))
    # else:
    #    plt.vlines(x=[start_time,delta_t+start_time], ymin=0, ymax=max(concentration_in))
    #    plt.annotate("$\Delta t$ = " + str(delta_t) + ", Average Outflow concentration has \nreached " + str(percentage) + " % of input concentraton.", (stop_time + 1 ,0.8*max(concentration_in)))

    ax2.legend(loc=0)
    plt.savefig(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.png",
        dpi=300,
    )
    np.savetxt(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.txt",
        avg_concentration,
    )
    #plt.show()

    return avg_concentration


def scenario_2b(task_id, task_root, single_file, var_name, normalize):
    """
    Dirac impulse of concentration for 1 day at a single point.
    """

    avg_concentration, times, time_steps, nodes = get_avg_conc_and_times_from_polyline(
        task_id=task_id, task_root=task_root, single_file=single_file
    )

    times_rfd, concentration_in = extract_rfd(task_root, rfd=1, export=True)

    start_time = get_start_time(np.asarray([times_rfd, concentration_in]).transpose())

    # index of maximum concentration - 2 (because of initial time step beeing 0)
    delta_t = times[np.argmax(avg_concentration)] - start_time
    stop_time = start_time + delta_t

    fig, ax1 = plt.subplots(figsize=(10, 7))
    ax1.set_title(os.path.basename(task_root))
    ax1.set_xlabel("days")
    ax1.set_ylabel("concentration in", color="b")
    ax1.tick_params("y", colors="b")
    ax1.plot(
        times / normalize,
        concentration_in,
        linewidth=4,
        label="input concentration",
        color="b",
    )
    ax1.legend(loc=2)
    plt.vlines(
        x=[start_time / normalize, stop_time / normalize],
        ymin=0,
        ymax=np.max(concentration_in),
        color="g",
    )
    plt.annotate(
        "  $\Delta t$ = " + str(delta_t / normalize),
        (stop_time / normalize, 0.3 * np.max(concentration_in)),
    )
    ax2 = ax1.twinx()
    ax2.plot(
        times / normalize,
        avg_concentration,
        linewidth=2,
        label="averaged outflow concentration",
        color="orange",
    )
    ax2.set_ylabel("concentration out", color="orange")
    ax2.tick_params("y", colors="orange")

    # if delta_t == " greater than modelling period...":
    #    plt.annotate(str(percentage) + " % is" + delta_t + " \nReached {:.3f} % in modelling period.".format(conc_out / conc_in), (np.median(times) + 1 ,0.8*max(concentration_in)))
    # else:
    #    plt.vlines(x=[start_time,delta_t+start_time], ymin=0, ymax=max(concentration_in))
    #    plt.annotate("$\Delta t$ = " + str(delta_t) + ", Average Outflow concentration has \nreached " + str(percentage) + " % of input concentraton.", (stop_time + 1 ,0.8*max(concentration_in)))

    ax2.legend(loc=0)
    plt.savefig(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.png",
        dpi=300,
    )
    np.savetxt(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.txt",
        avg_concentration,
    )
    #plt.show()


def scenario_2a(
    task_id, task_root, single_file, var_name, normalize, percentage
):
    """
    Script to plot the averaged concentration along a polyline vs the input concentration.
    """

    avg_concentration, times, time_steps, nodes = get_avg_conc_and_times_from_polyline(
        task_id=task_id, task_root=task_root, single_file=single_file
    )
    times_rfd, concentration_in = extract_rfd(task_root, rfd=1, export=True)
    start_time = get_start_time(np.asarray([times_rfd, concentration_in]).transpose())

    # 95 % of concentration
    for time, conc_out, conc_in in zip(times, avg_concentration, concentration_in):
        if conc_out >= percentage * conc_in and time > start_time:
            delta_t = time - start_time
            stop_time = delta_t + start_time
            print(
                "Avegerage concentration at the outlet has reached "
                + str(percentage)
                + " % of input concentration after "
                + str(delta_t / normalize)
                + " days."
            )
            break
        if time == np.max(times):
            print(
                "Concentration has not reached "
                + str(percentage)
                + " % of input concentration within modelling period. Reached {:.3f} %.".format(
                    conc_out / conc_in
                )
            )
            delta_t = " greater than modelling period..."

    plt.figure(figsize=(10, 7))
    plt.plot(
        times / normalize,
        avg_concentration,
        linewidth=4,
        label="averaged outflow concentration",
        color="orange",
    )
    plt.plot(
        times / normalize,
        concentration_in,
        linewidth=2,
        label="input concentration",
        color="blue",
    )
    plt.xlabel("days")
    plt.ylabel("concentration")
    plt.title(os.path.basename(task_root))
    if delta_t == " greater than modelling period...":
        plt.annotate(
            "  "
            + str(percentage)
            + " % is"
            + str(delta_t)
            + " \nReached {:.3f} % in modelling period.".format(conc_out / conc_in),
            (np.median(times) / normalize, 0.8 * max(concentration_in)),
        )
    else:
        plt.vlines(
            x=[start_time / normalize, delta_t / normalize + start_time / normalize],
            ymin=0,
            ymax=max(concentration_in),
            color="g",
        )
        plt.annotate(
            "  $\Delta t$ = "
            + str(delta_t / normalize)
            + ", Average Outflow concentration has \nreached "
            + str(percentage)
            + " % of input concentraton.",
            (stop_time / normalize, 0.8 * max(concentration_in)),
        )
    plt.legend()
    plt.savefig(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.png",
        dpi=300,
    )
    #plt.show()

    np.savetxt(
        task_root + "/" + os.path.basename(single_file) + "_avg_concentration.txt",
        avg_concentration,
    )

    return avg_concentration


if __name__ == "__main__":

    task_root = "/Users/houben/phd/modelling/20190602_Tc_vs_Tr/transport/setup/TEST2b_transport_stor_0.01_kf_0.001_disp_1_1_rech_1.16e-08_time_10000_864000_tort_0.5"
    task_id = "transect"
    parent_dir = "/Users/houben/phd/modelling/20190602_Tc_vs_Tr/transport/setup/third_tries/*/"
    file_name = "transect_ply_obs_01000_t5.tec"
    single_file = task_root + "/" + file_name
    var_name = "tracer"
    # normalize = 86400 if times are in seconds to end up with days
    normalize = 86400
    # - params for scenario_2a
    percentage = 0.95


    #scenario_2a(task_id, task_root, single_file, var_name, normalize, percentage)
    #scenario_2b(task_id, task_root, single_file, var_name, normalize)
    #scenario_2c(task_id,task_root, single_file, var_name, normalize, percentage, ref_conc=None)
    #scenario_2d(task_id, task_root, single_file, var_name, normalize)

    def run_all():
        import glob
        import os.path

        directories = glob.glob(parent_dir)
        for directory in directories:
            print("Working on dir: " + directory)
            folder = os.path.basename(directory[:-1])
            task_root = directory[:-1]
            single_file = task_root + "/" + file_name
            if folder[:2] == "2a":
                scenario_2a(
                    task_id, task_root, single_file, var_name, normalize, percentage
                )
            elif folder[:2] == "2b":
                scenario_2b(task_id, task_root, single_file, var_name, normalize)
            elif folder[:2] == "2c":
                scenario_2c(
                    task_id,
                    task_root,
                    single_file,
                    var_name,
                    normalize,
                    percentage,
                    ref_conc=None,
                )
            elif folder[:2] == "2d":
                scenario_2d(task_id, task_root, single_file, var_name, normalize)
            else:
                print("Nope, thats not an ogs model run folder: " + folder)
