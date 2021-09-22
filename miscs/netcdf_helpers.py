import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def extract_timeseries(path, xvalue, yvalue, varstring, xstring, ystring, timestring="time"):
    """
    Extracts and returns a time series based on coordinates XY from a netCDF file.

    Parameters
    ----------
    path : string or netCDF4.Dataset
        Path to netCDF file including file name.
        netCDF4.Dataset: if provided this dataset will 
    xvalue : float
        Value of the x-coordinate.
    yvalue : float
        Value of the y-coordinate.    
    varstring : string
        Name of the variable to extract.
    xstring : string
        Name of the x-coordinate.
    ystring : string
        Name of the y-coordinate.
    timestring : string, default: 'time'
        Name of the time variable

    Returns
    -------
    timeseries_df : dataframe
        Dataframe with two columns: ["Date", "timeseries"]
        Date: Date time object
        timeseries: time series of the variable

    """
    netcdffile = netCDF4.Dataset(path, "r")
    rechts_arr = netcdffile.variables[xstring][:]
    hoch_arr = netcdffile.variables[ystring][:]
    times = netcdffile.variables[timestring]
    # convert numerics to date
    jd = netCDF4.num2date(times[:], times.units)

    # function to find index to nearest point
    def near(array, value):
        diff_array = np.abs(array - value)
        id = np.argmin(diff_array)
        # id = np.unravel_index(np.argmin(diff_array), diff_array.shape)
        return id
    # find nearest point to desired location
    ix = near(rechts_arr, xvalue)
    iy = near(hoch_arr, yvalue)
    # find value of projected coordinate system of nearest location
    rechts_out = netcdffile.variables[xstring][ix]
    hoch_out = netcdffile.variables[ystring][iy]
    # create dataframe from data to save it as txt (can be improved)
    timeseries = netcdffile.variables[varstring][:, iy, ix]
    time_df = pd.DataFrame(jd, columns=["Date"])
    time_df["Date"] = pd.to_datetime(time_df["Date"].dt.date, format="%Y-%m-%d")
    time_df["Date"] = time_df["Date"].dt.date
    timeseries_df = pd.DataFrame((timeseries), columns=["timeseries"])
    time_timeseries_df = pd.concat([time_df, timeseries_df], axis=1)
    time_timeseries_df.set_index("Date", inplace=True)
    # time_timeseries_df.to_csv(
    #     directory + "/" + str(rechts_out) + "_" + str(hoch_out) + "_" + name + "_timeseries.txt",
    #     header=None,
    #     sep=" ",
    #     index=None,
    # )
    return time_timeseries_df


def extract_and_average(path, xvalues, yvalues, varstring, xstring, ystring, timestring="time"):
    """
    Extract time series for a list of xy values and average the resulting series.
    
    Parameters
    ----------
    path : string
        Path to netCDF file including file name.
    xvalues : list of floats
        Value of the x-coordinate.
    yvalues : list of floats
        Value of the y-coordinate.    
    varstring : string
        Name of the variable to extract.
    xstring : string
        Name of the x-coordinate.
    ystring : string
        Name of the y-coordinate.
    timestring : string, default: 'time'
        Name of the time variable

    Returns
    -------
    date : list
        List of datetime objects.
    timeseries : list
        List of extracted time series.

    """
    
    timeseries_df_list = []
    i = 0
    for i, (xvalue, yvalue) in enumerate(zip(xvalues, yvalues)):
        timeseries_df = extract_timeseries(path, xvalue, yvalue, varstring, xstring, ystring, timestring="time")
        timeseries_df.rename(columns={"timeseries":(str(xvalue)+"_"+str(yvalue))}, inplace=True)
        timeseries_df_list.append(timeseries_df)
        # if i == 3:
        #     break
        # i += 1
        print(str(i) + "/" + str(len(xvalues)))
    average_df = pd.concat(timeseries_df_list, axis=1)

    print("Ready")
    
    return average_df




if __name__ == "__main__":
    pass


    # Testing the function extract_timeseries
    # ---------------------------------------
    # path = "/Users/houben/phd/mHM/GW_recharge_main/raw_data/mHM_recharge_Main_basin_1955-2019_daily.nc"
    # varstring = "recharge"
    # xvalue = 4387322.53774309
    # yvalue = 5517587.77908758
    # xstring="easting"
    # ystring="northing"
    # date, timeseries = extract_timeseries(path, xvalue, yvalue, varstring, xstring, ystring)


    # # Testing extract_and_average
    # # ---------------------------
    # coordinates_df = pd.read_csv("/Users/houben/phd/studies/SA-Main-Toy/gis/slices/slice_2_xy.csv")
    # path = "/Users/houben/phd/mHM/GW_recharge_main/raw_data/mHM_recharge_Main_basin_1955-2019_daily.nc"
    # xvalues = coordinates_df.X.tolist()
    # yvalues = coordinates_df.Y.tolist()
    # varstring = "recharge"
    # xstring="easting"
    # ystring="northing"
    # average_df = extract_and_average(path, xvalues, yvalues, varstring, xstring, ystring, timestring="time")
    # average_df.to_csv("/Users/houben/phd/studies/SA-Main-Toy/recharge/slice_2_xy_recharge_df.csv")
    # average_df = pd.read_csv("/Users/houben/phd/studies/SA-Main-Toy/recharge/slice_2_xy_recharge_df.csv")
        
    # get the mean for each row
    slice_2_mean_recharge = average_df.mean(axis=1)
