# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
# ------------------------
# set matplotlib style
CWD = os.getcwd()
CWD_p = os.path.dirname(CWD)
plt.style.use(["ggplot", os.path.join(CWD_p, "miscs", "mpl_style_1.mplstyle")])

def plot_wn_mhm_recharge(wn, mhm, times, savepath):
    factor = 1000 * 86400
    
    wn = wn * factor
    mhm = mhm * factor

    times = times / 86400 / 365
    plt.figure(figsize=(16,9))
    plt.plot(times, wn, label="White Noise-Like", linewidth=0.1)
    plt.plot(times, mhm, label="mHM", linewidth=1)

    plt.legend()
    plt.ylabel("Recharge [$mmd^{-1}$]")
    plt.xlabel("Time [years]")
    
    plt.savefig(savepath + "/recharge.pdf", dpi=300)
    



if __name__ == "__main__":

    # plot recharge for paper spectral analysis
    times = np.loadtxt("/Users/houben/phd/paper/spectral_anaysis_1/recharge/rfd_curve#1_x_values.txt")
    mhm = np.loadtxt("/Users/houben/phd/paper/spectral_anaysis_1/recharge/rfd_curve#1_y_values.txt")
    wn = np.loadtxt("/Users/houben/phd/paper/spectral_anaysis_1/recharge/20190304_recharge_daily.txt", usecols=1)
    savepath = "/Users/houben/phd/paper/spectral_anaysis_1/recharge"
    plot_wn_mhm_recharge(wn, mhm, times, savepath)

