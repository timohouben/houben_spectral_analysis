import matplotlib.pyplot as plt
import numpy as np
from ogs5py.reader import readtec_polyline
import sys
sys.path.append("/Users/houben/phd/python/scripts/spectral_analysis")
from ogs_helpers.transect_plot import extract_timeseries
import os
"""
Script to check if the extract time series works.
Problem: txt files which were exported with keyword "mean" were wrong,
at least partly...
"""
# project directory
path = "/Users/houben/phd/studies/spectral_analysis/20210823_spectral_analysis_homogeneous_neu_new_mesh/2079_769.7771_3.98e-03_2.00e-05"

# ----- extract a single time series manually
# task_id = "transect"
# task_root = path
# process = "GROUNDWATER_FLOW"
# obs = "obs_00500"
# tecs = readtec_polyline(task_id=task_id, task_root=path)
# time_steps = len(tecs["GROUNDWATER_FLOW"][obs]["TIME"])
# number_of_columns = tecs[process][obs]["HEAD"].shape[1]

# time = tecs["GROUNDWATER_FLOW"]["obs_00500"]["TIME"]


# head_ogs_timeseries_each_obs = []
# for step in range(10951):
#     # calculates the mean of each time step
#     head_ogs_timeseries_each_obs.append(
#         np.mean(tecs[process][obs]["HEAD"][step, :])
#     )

# plt.plot(head_ogs_timeseries_each_obs)
# plt.show()



# ----- extract the time series with the function

extract_timeseries(path, which="max", process="GROUNDWATER_FLOW")

# ---- Plot the extracted time series
min500 = np.loadtxt(os.path.join(path, "head_ogs_obs_00600_min.txt"))
max500 = np.loadtxt(os.path.join(path, "head_ogs_obs_00600_max.txt"))
mean500 = np.loadtxt(os.path.join(path, "head_ogs_obs_00600_mean.txt"))

plt.plot(min500, label="min")
plt.plot(max500, label="max")
plt.plot(mean500, label="mean")
plt.legend()
plt.show()