# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division

# ------------------------------------------------------------------------------

def get_kf_from_blocks(path):
    """
    This function returns a list with kf values from the mmp blocks in the
    ogs .mmp-file. It's sorted like the blocks are sorted, i.e. like the material
    groups.
    """
    from ogs5py import OGS
    ogsmodel = OGS(task_root=path)
    ogsmodel.load_model(task_root=path)
    kf_list = []
    number_of_blocks = ogsmodel.mmp.get_block_no()
    for i in range(number_of_blocks):
        kf_list.append(ogsmodel.mmp.get_block(i)['PERMEABILITY_TENSOR'][0][1])
    return kf_list


def get_ogs_parameters(path, noKf=False):
    """Return Ss, kf from the first block of the .mmp-file
    Return time_step_size and time_steps from the .tim-file."""


    from ogs5py import OGS

    ogsmodel = OGS(task_root=path)
    ogsmodel.load_model(task_root=path)
    # dir(ogsmodel.mmp)
    # print(vars(ogsmodel.mmp))
    # for any reason is this configuration not working with python2
    # Ss = ogsmodel.mmp.get_block()['STORAGE'][0][1]
    # kf = ogsmodel.mmp.get_block()['PERMEABILITY_TENSOR'][0][1]
    # print(vars(ogsmodel.mmp))

    # this configuration is probably only working for this kind of OGS setups
    # because it is indexing the contents of the mmp file and not refering to a
    # dictionary with keys!!!
    # alternative: mpd.get_block(0)[""]
    Ss = float(ogsmodel.mmp.cont[0][1][0][1])
    kf = float(ogsmodel.mmp.cont[0][2][0][1])
    time_step_size = float(ogsmodel.tim.cont[0][3][0][1])
    time_steps = float(ogsmodel.tim.cont[0][3][0][0])
    if noKf == True:
        return Ss, time_step_size, time_steps
    else:
        return Ss, kf, time_step_size, time_steps


if __name__ == "__main__":
    #path = "/Users/houben/PhD/modelling/20190304_spectral_analysis_homogeneous/models/100_sample2_351_1.10e-05_1.00e-03"
    #print(get_ogs_parameters(path))
    print(get_kf_from_blocks("/Users/houben/Desktop/eve_work/20190917_generate_ogs_layered_ensemble/setup/1001_reali_0_stor_0.01_rech_whitenoise", 30))
