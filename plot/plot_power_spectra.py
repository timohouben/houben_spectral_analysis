# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------


def plot_spectrum(
    data,
    frequency,
    name=None,
    labels=None,
    path=None,
    lims=None,
    linestyle="-",
    marker="",
    #markersize=None,
    grid="both",
    unit="[Hz]",
    heading="None",
    figtxt=None,
    comment="",
    r2=None,
    ext=".png",
    **kwargs
):
    """
    Function to plot one or multiple power spectra.

    Parameters
    ----------
    data : 2-D array
        Each row represents a seperate power spectrum.
    frequency : 1-D array
        Corresponding frequencies of data.
    name : string
        Name of file. If None, time is used.
    labels : X item list
        Labels for different power spectra as list in same order as data.
    path : string
        Path to store the image.
    lims : list with 2 tuples
        lims[0] = x limit as tuple (xmin,xmax)
        lims[1] = y limit as tuple (ymin,ymax)
        e.g. lims = [(1e-8,1e-4),(1e0,1e5)]
    linestyle : X item list
        List with linestyles for differenct spectra.
    marker : X item list
        List with marker for differenct spectra.
    grid : string
        "major", "minor", "both", "none"
    unit : string
        Unit of frequency.
    heading : string
        Provide a heading for the image. If None, no heading.
    figtxt : string (multiline possible)
        Provide an annotation for a box below the figure. If None, no annotaion.
    comment : string
        A comment.
    r2 : float
        An r2 value from the fit analytical sol. to data.
    ext : string
        File extension of the saved image.
        '.png', '.eps', '.jpg'

    Yields
    ------
    One saved image in path.
    """

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.style.use("ggplot")
    font = {"family": "DejaVu Serif", "weight": "normal", "size": 20}
    plt.rc("font", **font)
    plt.rc("legend", fontsize=20)
    plt.grid(which=grid, color="grey", linestyle="-", linewidth=0.2)
    plt.figure(figsize=[20, 10], dpi=300)
    if np.ndim(data) == 1:
        plt.loglog(
            frequency,
            data,
            label=str(labels[0]),
            linewidth=1,
            linestyle=linestyle,
            marker=marker,
            **kwargs
            #markersize=markersize,
        )
    elif (np.ndim(data) != 1) & (np.ndim(frequency) == 1):
        for i, spectrum in enumerate(data):
            plt.loglog(
                frequency,
                spectrum,
                label=labels[i],
                linewidth=1,
                linestyle=linestyle[i],
                marker=marker[i],
                #markersize=markersize[i],
                **kwargs
            )
    else:
        for i, spectrum in enumerate(data):
            plt.loglog(
                frequency[i],
                spectrum,
                label=labels[i],
                linewidth=1,
                linestyle=linestyle[i],
                marker=marker[i],
                #markersize=markersize[i],
                **kwargs
            )
    if lims is not None:
        plt.xlim(lims[0])
        plt.ylim(lims[1])
    if heading is not None:
        plt.title(heading)
    if labels is not None:
        plt.ylabel("Power Spectrum")
        plt.xlabel("Frequency %s" % unit)
    plt.legend(loc="best")
    if figtxt is not None:
        plt.figtext(
            0.135,
            -0.1,
            figtxt,
            horizontalalignment="left",
            bbox=dict(boxstyle="square", facecolor="#F2F3F4", ec="1", pad=0.4, alpha=1),
        )
    if r2 is not None:
        if lims is None:
            plt.text(np.min(frequency)*10, np.min(data)*10, r'$R^2\ =\ {:.4f}$'.format(r2))
        else:
            plt.text(np.min(lims[0][0])*10, np.min(lims[1][0])*10, r'$R^2\ =\ {:.4f}$'.format(r2))
    if path is not None:
        import datetime

        if name is None:
            name = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        plt.savefig(
            path + "/" + comment + name + ext, pad_inches=0.5, bbox_inches="tight"
        )
    else:
        plt.show()
    plt.close()


def plot_shh_anal_loc():
    """
    Function to plot multiple analytical power spectra along e.g. an aquifer in a 3D plot.
    Still not working because plot3d has issues with log scale ...
    """

    # set parameters
    data_points = 8000
    time_step_size = 86400
    aquifer_length = 1000
    Sy = 1e-1
    T = 0.001
    from calc_tc import calc_tc

    tc = calc_tc(aquifer_length, Sy, T)
    print(tc)

    import sys

    # add search path for own modules
    sys.path.append("/Users/houben/PhD/python/scripts/spectral_analysis")
    from shh_analytical import shh_analytical
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D
    import scipy.fftpack as fftpack

    # create an input signal
    np.random.seed(123456789)
    input = np.random.rand(data_points)
    spectrum = fftpack.fft(input)
    # erase first half of spectrum
    spectrum = abs(spectrum[: round(len(spectrum) / 2)]) ** 2
    spectrum = spectrum  # [1:]
    # X contains the different locations
    X = np.linspace(0, aquifer_length - 1, 10)
    X = [100, 900]
    # Y contains the frequencies, erase first data point because it's 0
    Y = abs(fftpack.fftfreq(len(input), time_step_size))[: round(len(input) / 2)]
    Y = np.log10(Y[1:])
    Z = np.zeros((len(Y), len(X)))
    for i, loc in enumerate(X):
        Z[:, i] = np.log10(
            shh_analytical(
                (Y, spectrum), Sy=Sy, T=T, x=loc, L=aquifer_length, m=5, n=5, norm=False
            )
        )
    # erase first data point from Z for each location
    print(Z)
    print(Y)
    import matplotlib.pyplot as plt

    plt.plot(Y, Z[:, 0])
    # plt.loglog(Y,Z[:,0])
    #    X, Y = np.meshgrid(X, Y)
    #    fig = plt.figure()
    #    ax = Axes3D(fig)
    # surf = ax.plot_surface(
    #    X, Y, Z, rstride=1, cstride=2, shade=False, linewidth=1, cmap="Spectral_r"
    # )
    #    surf = ax.plot_wireframe(X, Y, Z, rstride=0, cstride=1, cmap=cm.magma)
    # surf.set_edgecolors(surf.to_rgba(surf._A))
    # surf.set_facecolors("white")
    # ax1 = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=0)

    #    ax.set_xlabel("Location [m]")
    #    ax.set_ylabel("Frequency [Hz]")
    #    ax.set_zlabel("log Spectral Density")

    # ax.set_zscale("log")
    # ax.yaxis.set_scale("log")
    # ax.zaxis._set_scale('log')
    # ax.set_yscale("log")
    plt.show()


if __name__ == "__main__":
    plot_shh_anal_loc()
    # plot_shh_anal_S(aquifer_length=1000, time_step_size=86400)

    # Test for function plot_spectrum
#    from power_spectrum import power_spectrum
#    import numpy as np

#    frequency, data1 = power_spectrum(
#        np.random.rand(1000), np.random.rand(1000), 86400, o_i="o"
#    )
#    frequency, data2 = power_spectrum(
#        np.random.rand(1000), np.random.rand(1000), 86400, o_i="o"
#    )
#    frequency, data3 = power_spectrum(
#        np.random.rand(1000), np.random.rand(1000), 86400, o_i="o"
#    )
# data = np.vstack((data1, data2, data3))
# labels = ["head1", "head2", "head3"]
# linestyle = ["-", "--", ":"]
# path = "/Users/houben/Desktop/TEST"
# plot_spectrum(
#        data,
#        frequency,
#        labels,
#        path,
#        lims=None,
#        figtxt="nonoasndoand\nasdasd",
#        linestyle=linestyle,
#    )
