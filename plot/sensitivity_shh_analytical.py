# %matplotlib qt
import sys

sys.path.append("/Users/houben/phd/python/scripts/spectral_analysis")
from analysis.shh_analytical import shh_analytical
import numpy as np
import matplotlib.pyplot as plt
from analysis.power_spectrum import power_spectrum
from analysis.calc_tc import calc_tc
import os
import string
"""
Script to plot the shh_analytical for different parameters for paper
spectral analysis 1.
"""
# ------------------------
# matplotlib costum design
# ------------------------
# set matplotlib style
CWD = os.getcwd()
CWD_p = os.path.dirname(CWD)
plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])
# plt.rc('text', usetex=True)
#plt.rc("font", family="serif")
#plt.rc("xtick", labelsize="x-large")
#plt.rc("ytick", labelsize="x-large")
# fig = plt.figure(figsize=(16, 9))
# ------------------------
# path to save the image
savepath = "/Users/houben/phd/paper/spectral_anaysis_1/plots"

# time step size in seconds
time_step_size = 86400
# load recharge
recharge = np.loadtxt(
    "/Users/houben/phd/paper/spectral_anaysis_1/recharge/20190304_recharge_daily.txt"
)
# get the spectrum of the recharge
frequency, sww = power_spectrum(recharge[:, 1], recharge[:, 1], time_step_size, o_i="i")
# shh_analytical(X, Sy, T, x, L, m=None, n=None, norm=False, convergence=0.01)
X = (frequency, sww)
# list for Storativities
S_list = [0.3, 0.2, 0.1, 0.01, 0.001]
# list for Transmissivitites
T_list = [0.01, 0.001, 0.0001, 0.00001, 0.000001]
# list for locations
x_list = [10, 100, 400, 600, 900]
# list for lengths
L_list = [1000, 2000, 5000, 10000, 20000]
# define standard values for remaining parameters
S_std, T_std, x_std, L_std = 0.01, 0.001, 500, 1000
# rows, cols for axis
rows, cols = 2, 2
# figsize
figsize = (16, 9)
# linestyle
ls = "-"
# linewidth
lw = "0.5"
# markerstyle
marker = " "
# markersize
ms = "2"
# add vertical lines for 1 year, 1 month and 1 week
v_lines = [1/(365*86400), 1/(30*86400), 1/(7*86400)]
# define labels for v lines
v_line_labels = ["Year", "Month", "Week"]
# define a list for numbers for plots
alpha_list = ["({})".format(i) for i in string.ascii_lowercase] + [
    "({})".format(i) for i in string.ascii_uppercase
]


axs = plt.figure(figsize=figsize, constrained_layout=True).subplots(
    rows, cols, sharex=True, sharey=True
)
for S in S_list:
    tc = calc_tc(L_std, S, T_std, which="dupuit")
    axs[0, 0].loglog(
        frequency,
        shh_analytical(X, T_std, S, x_std, L_std),
        label=fr"S = {S:.0e}, $t_c$ = {tc:.0f} d", linestyle=ls, linewidth=lw, marker=marker, ms=ms
    )
    axs[0, 0].legend()
for T in reversed(T_list):
    tc = calc_tc(L_std, S_std, T, which="dupuit")
    axs[0, 1].loglog(
        frequency,
        shh_analytical(X, T, S_std, x_std, L_std),
        label=fr"T = {T:.0e} $m^2s^{-1}$, $t_c$ = {tc:.0f} d", linestyle=ls, linewidth=lw, marker=marker, ms=ms
    )
    axs[0, 1].legend()
for x in x_list:
    tc = calc_tc(L_std, S_std, T_std, which="dupuit")
    axs[1, 0].loglog(
        frequency,
        shh_analytical(X, T_std, S_std, x, L_std),
        label=fr"x = {x} m, $t_c$ = {tc:.0f} d", linestyle=ls, linewidth=lw, marker=marker, ms=ms
    )
    axs[1, 0].legend()
for L in reversed(L_list):
    tc = calc_tc(L, S_std, T_std, which="dupuit")
    axs[1, 1].loglog(
        frequency,
        shh_analytical(X, T_std, S_std, x_std, L),
        label=fr"L = {L} m, $t_c$ = {tc:.0f} d", linestyle=ls, linewidth=lw, marker=marker, ms=ms
    )
    axs[1, 1].legend()
axs[0, 0].set_ylabel("Power Spectrum")
axs[1, 0].set_ylabel("Power Spectrum")
axs[1, 0].set_xlabel("Frequency [Hz]")
axs[1, 1].set_xlabel("Frequency [Hz]")

for i, ax in enumerate(axs.flat):
    for v_line, v_line_label in zip(v_lines, v_line_labels):
        ax.axvline(x=v_line, linestyle="--", color="grey", alpha=0.2)
        ax.annotate(
                    v_line_label, (v_line * 1.1, 1e-5), alpha=0.5, fontsize="16", color="grey"
                )
    # make an annotation for numbering of subplots
    ax.annotate(alpha_list[i], (1e-7, 1e7), fontsize=16)
plt.savefig(os.path.join(savepath, "sensitivity_shh_analytical.pdf"))