# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division
import numpy as np
# ------------------------------------------------------------------------------

def aggregate_series(series, n):
    """
    Takes the mean of n entries in a series. If len(series) % n != 0 remaining
    values will be popped.

    Parameters
    ----------

    series : 1D array
        Series to aggregate.
    n : integer
        Number of entries to aggregate to.
    """
    import numpy as np

    series_agg = []
    len_series_agg = len(series) // n
    for i in range(len_series_agg):
        mean_n = np.mean(series[n*i:n*(i+1)])
        series_agg.append(mean_n)
        #print(mean_n)
    series_agg = np.asarray(series_agg)
    print("You aggregated a time series with n = ", str(n))
    return series_agg

def identify_numbers_from_string(string, index=None):
    """
    Pass a string or a list of strings and return a list of numbers in this string.

    Parameters
    ----------
    string : string or list
        Give the string which should be searched for numbers of any kind.
        Lists will be concatinated.
    index : integer
        Return only one value with index.

    """
    import re
    if not isinstance(string, str):
        string = np.ravel(string)
        string = ",".join([str(i) for i in string])

    if index is None:
        found = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", string)
        if found == None:
            return np.nan
        else:
            return found
    elif index is not None:
        found = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", string)[index]
        if found == None:
            return np.nan
        else:
            return found


def detrend(timeseries, type="linear"):
    """
    Function to detrend the time series.

    Parameters
    ----------

    timeseries : 1D array
        time series to be detrendet
    type : string
        linear or ...?

    Yields
    ------

    timeseries : 1D array
        Detrendet time series.
    """
    from scipy.signal import detrend

    return detrend(timeseries, type=type)


def cut(x, y, threshold, keep="before"):
    """
    Cut a series x at given threshold and cut corresponding series y
    at same index. Threshold will be kept in series when 'before' is chosen.
    If 'after' was chosen, threshold will be not in list any more.

    Parameters
    ----------
    x : 1D array
        Series to which the threshold refers. Should be sorted either
        ascending or descending.
    y : 1D array
        Series which will be cut at same index as x.
    threshold : float
        Threshold where x series will be cut (excluding)
    keep : string
        'before' : values with index smaller than the threshold will be kept
        'after' : values with index after threshold will be kept

    Yields
    ------
    x_cut : 1D array, list
        Cut series.
    y_cut : 1D array, list
        Cut series.
    """

    import numpy as np

    if np.shape(x) != np.shape(y):
        raise ValueError
        print("x and y must have same length.")
    if np.asarray(x).ndim != 1:
        raise ValueError
        print("x and y must have dimension = 1.")

    if [i for i in sorted(x)] == [i for i in x]:
        if threshold < x[0]:
            raise ValueError
            print("Your threshold is to low. Not cutting list.")
        if threshold > x[-1]:
            raise ValueError
            print("Your threshold is to high. Not cutting list.")
        for i, item in enumerate(x):
            if item > threshold:
                if keep == "before":
                    return x[:i], y[:i]
                elif keep == "after":
                    return x[i:], y[i:]
    elif [i for i in sorted(x, reverse=True)] == [i for i in x]:
        if threshold > x[0]:
            raise ValueError
            print("Your threshold is to high. Not cutting list.")
        if threshold < x[-1]:
            raise ValueError
            print("Your threshold is to low. Not cutting list.")
        for i, item in enumerate(x):
            if item < threshold:
                if keep == "before":
                    return x[:i], y[:i]
                elif keep == "after":
                    return x[i:], y[i:]
    else:
        raise ValueError(
            "Your series x is not sorted. Sort it either ascending or descending."
        )


def percent_fraction(a, b):
    """
    Returns fraction of a in b as %.
    a / b * 100

    Parameters
    ----------
    a : float

    b : float

    """
    return a / b * 100


def percent_difference_fraction(a, b):
    """
    Returns fraction of difference of a and b to a as %.
    (a - b)/a * 100

    Parameters
    ----------
    a : float

    b : float

    """
    return (a - b) / a * 100


def percent_difference_fraction_log(a, b):
    """
    Returns fraction of difference of log10(a) and log10(b) to log10(a) as %.
    ( log10(a) - log10(b) ) / log10(a) * 100

    Parameters
    ----------
    a : float

    b : float

    """
    import numpy as np

    return (np.log10(a) - np.log10(b)) / np.log10(a) * 100


def difference_fraction_log(a, b):
    """
    Returns fraction of difference of log10(a) and log10(b) to log10(a).
    ( log10(a) - log10(b) ) / log10(a) * 100

    Parameters
    ----------
    a : float

    b : float

    """
    import numpy as np

    return (np.log10(a) - np.log10(b)) / np.log10(a)


def combine_results(path_to_multiple_projects, name):
    """
    Searches for results.csv in all subdirectories and combines these files. All csv files must have the same header!

    Parameters
    ----------

    path_to_multiple_projects : string
        Path in which you want to start searching. Results will be stored on this level.
    name : string
        If file name contains this sting it will be selected for combination.
    """

    import os
    import os.path

    if os.path.exists(path_to_multiple_projects + "/combined_results") == False:
        os.mkdir(path_to_multiple_projects + "/combined_results")

    file_paths = []
    for dirpath, dirnames, filenames in os.walk(path_to_multiple_projects):
        for filename in [f for f in filenames if name in f]: 
            file_paths.append(dirpath + "/" + str(filename))

    # get header from first file in list
    with open(file_paths[0]) as f:
        header = f.readline()
        f.close()

    csv_merge = open(
        path_to_multiple_projects + "/combined_results" + "/" + filename[:-4] + "_merge.csv", "w"
    )
    csv_merge.write(header)

    for file in file_paths:
        csv_in = open(file)
        for line in csv_in:
            if line.startswith(header):
                continue
            csv_merge.write(line)
        csv_in.close()
    csv_merge.close()
    print("Created consolidated CSV file : " + name[:-4] + "_merge.csv")


def get_filename_from_rfd_top_com(
    path="/Users/houben/phd/modelling/20190430_spectral_analysis_heterogeneous/models/var_1_len_100_mean_-10_seed_12345_results",
    t_id="transect",
):
    rfd_file = open(path + "/" + t_id + ".rfd", "r")
    for line in rfd_file:
        position_slash = line.rfind("/")
        position_newline = line.rfind("\n")
        break
    return line[position_slash + 1 : position_newline]





if __name__ == "__main__":
    pass
    # test for combine_results()
    # combine_results("/Users/houben/Desktop/TEST_combine")

    # # test for def aggregate_series(series, n)
    # # ----------------------------------------
    # import numpy as np
    # np.random.seed(199)
    # series = np.random.randint(0,100,100)
    # n = 99
    # test_array = []
    # for i in range(len(series) // n):
    #     test_array.append(np.mean(series[i*n:(i+1)*n]))
    # test_array = np.asarray(test_array)
    # series_agg = aggregate_series(series, n)
    # if np.all(test_array == series_agg):
    #     print("Test bestanden!")
    # else:
    #     print("Test nicht bestanden!")

    # # Test identify_numbers_from_string
    string = ["s15454353", "as1d"]
    string = [["s15454353"], "as1d"]
    string = np.array([["inf", "inf"],["inf", "inf"]])
    print(identify_numbers_from_string(string))
