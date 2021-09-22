# script to check the magnitude of the recharge spectra

# import modules
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# add search path for own modules
sys.path.append("/Users/houben/PhD/python/scripts/spectral_analysis")
# add search path for own modules on eve
sys.path.append("/home/houben/python/scripts/spectral_analysis")
# own modules
from analysis.power_spectrum import power_spectrum
from plot.plot_power_spectra import plot_spectrum
from analysis.shh_analytical import shh_analytical_fit, shh_analytical
from analysis.calc_tc import calc_tc
from analysis.processing import detrend

recharge_df = pd.read_csv(
            "/Users/houben/phd/studies/application_spectral_analysis/main/20200710_spectral_analysis/input/recharge/4416000.0_5522000.0_birkach_recharge_2019.txt",
            sep=" ",
            header=None,
            names=["date", "head"],
        )
recharge = recharge_df["head"].tolist()
time_step_size = 86400

frequency_output, Shh = power_spectrum(
    input=recharge,
    output=recharge,
    time_step_size=time_step_size,
    method="scipyffthalf",
    o_i="o",
)


plt.loglog(frequency_output, Shh)
plt.show()