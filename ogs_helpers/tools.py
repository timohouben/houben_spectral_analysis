def get_ogs_polyline_length(path, name):
    """
    Returns the length of a polyline. The number of points which make up the
    polyline must not be more than 2!

    Parameters
    ----------
    path : string
        Path to OGS project folder.
    name : string
        Name of the polyline which length should be returned.
    """

    from ogs5py import OGS
    from scipy.spatial.distance import euclidean
    found = False
    ogsmodel = OGS(task_root=path)
    ogsmodel.load_model(task_root=path)
    for line in ogsmodel.gli.POLYLINES:
        if line["NAME"] == name:
            found = True
            points = line["POINTS"]
            if len(points) > 2:
                print("The polyline has more than 2 points!")
                print("Exit!")
                break
            return euclidean(ogsmodel.gli.POINTS[points[0]],ogsmodel.gli.POINTS[points[1]])
        else:
            pass

    if found == False:
        raise NameError("No polyline with the name " + name + " has been found!")
    else:
        pass

def get_ogs_folders(path):
    """
    Returns a list of directories where OGS model runs have been setup based on the following file types. It does not decide whether the model has run or not.

    file_extensions_list = ["*.gli", "*.msh", "*.out", "*.pcs", "*.num"]

    Parameters
    ----------
    path : string
        Path to multiple OGS project folders.

    Yields
    ------
    project_folder_list : list of strings
        Containing all folder names where OGS runs have been set up.
    """

    import os
    import glob

    file_extensions_list = ["*.gli", "*.msh", "*.out", "*.pcs", "*.num"]
    project_folder_list = [f for f in os.listdir(str(path)) if not f.startswith(".")]
    project_folder_list_temp = project_folder_list.copy()

#    print(project_folder_list)
    for folder in project_folder_list:
#        print(folder)
        check_extensions = []
#        print(len(project_folder_list))
        for extension in file_extensions_list:
            if glob.glob(path + "/" + folder + "/" + extension):
                check_extensions.append(True)
#                print(extension + " in " + folder)
            else:
                check_extensions.append(False)
#                print(extension + " not in " + folder)
        if any(check_extensions) == False:
            project_folder_list_temp.remove(folder)
    return project_folder_list_temp


def get_ogs_task_id(path):
    """
    Grabs the name of the ogs task id from .bc file.

    Parameters
    ----------

    path : strig
        Path to ogs directoy.

    """
    import glob
    # glob the name of the ogs run
    string = str(glob.glob(path + "/*.bc"))
    pos1 = string.rfind("/")
    task_id = string[pos1 + 1 : -5]

    return task_id


if __name__ == "__main__":
    pass
    #name = "right"
    #path = "/Users/houben/phd/modelling/20190712_Constant_vs_Constant_Neumann/1b_steady_gw_flow_stor_0.1_kf_0.0001_C_1"
    #test = get_ogs_polyline_length(path, name)
