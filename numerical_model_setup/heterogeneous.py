#!/usr/bin/env python3
# -*- coding: utf-8 -*

## This script generates OGS setups for a heterogeneous, 2D, vertical model domain. Every model setup ist first generated with a steady state configuration, afterwards the OGS run starts, output is redirected into the "output folder" and input files are changed to input files with transient settings. THIS IS NOT YET WORKING FINE!!! YOU NEED TO ADJUST THE PCS FILE MANUALLY (ERASE THE KEYWORD STEADY).
# Depending on where you start the script (local, EVE) you have to set the path to the OGS executable first and also the CWD and the path to the recharge files!!!


import sys
import os
import numpy as np
from ogs5py import OGS, MPD
import shutil
from ogs5py.reader import readpvd, readtec_point
from gstools import SRF, Gaussian
import matplotlib.pyplot as plt

# get the current working directory
# CWD = os.getcwd()
CWD = "/work/houben/20190513_spec_anal_hetero_ensemble_1/"
# the name of this script
file_name = "20190231_generate_ogs_heterogeneous_ensemble.py"
# ------------------------domain configuration-------------------------------- #
length = 1000
thickness = 30
n_cellsx = int(length)
n_cellsz = int(thickness)
# ------------------------time configuration---------------------------------- #
time_start = 0
time_steps = np.array([365*30])
step_size = np.array([86400])
time_end = np.sum(time_steps*step_size)
# ------------------------heterogeneous field configuration------------------- #
dim=2
var_list=[1]#,5,10]
len_scale_list=[5, 15, 100, 500]
mean_list=[-10]#, -10, -12]
# This list was generated with np.random.randint(1000,10000,200)
seed_list=[9631, 2723, 6039, 6003, 7743, 8375, 8031, 6689, 3521, 6576, 8077, 6399, 5547, 3342, 1500, 5133, 1823, 2619, 4710, 5321, 4816, 2949, 2409, 9705, 9885, 5975, 9522, 2119, 9058, 3103, 9813, 6466, 7012, 1547, 5864, 9895, 8987, 8382, 2239, 5745, 7391, 2096, 5595, 8846, 2764, 6969, 9261, 2440, 5948, 3254, 1538, 2956, 7594, 4096, 2419, 8450, 5643, 9462, 7218, 2362, 5517, 7654, 8690, 6746, 1877, 6462, 6824, 2854, 5381, 7227, 5284, 9396, 5327, 4112, 5751, 7209, 2465, 3104, 3850, 3736, 6767, 4357, 4610, 3440, 9542, 2334, 2908, 5315, 5382, 5249, 6930, 4648, 1522, 6340, 5128, 8688, 5206, 6622, 6805, 4093, 1329, 8158, 4643, 1774, 9727, 6386, 9952, 8210, 9256, 3512, 6808, 6309, 6339, 7437, 6937, 8298, 6215, 7842, 2472, 8016, 7555, 4760, 7516, 8927, 4166, 2204, 9107, 2417, 3469, 8314, 7560, 8492, 9733, 2950, 6888, 4518, 8624, 2645, 8376, 2114, 3849, 3466, 6750, 1448, 5458, 4256, 2360, 7428, 9115, 6451, 1680, 7980, 5562, 4983, 3960, 2202, 9504, 8606, 9182, 7984, 2802, 3899, 5267, 5481, 4521, 8729, 5942, 3746, 4222, 8784, 9035, 2221, 6544, 5842, 4206, 3851, 5160, 4934, 8365, 3167, 2755, 4368, 9756, 2617, 5443, 1956, 7227, 4404, 7574, 2594, 5361, 4165, 5173, 9628, 3346, 5297, 7882, 4820,  6515, 7312]
# ------------------------ogs configuration---------------------------------- #
# name of the folder of one single ogs run
storage_list = [0.001]#, 0.1]
recharge_list = ['/work/houben/20190513_spec_anal_hetero_ensemble_1/recharge_daily.txt']#, '/work/houben/20190513_spec_anal_hetero_ensemble_1/recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt']
rech_abv_list = ['whitenoise']#, 'mHM']
# index for model runs
overall_count = 1000

for seed in seed_list:
    for var in var_list:
        for len_scale in len_scale_list:
            for mean in mean_list:
                for storage in storage_list:
                    for recharge, rech_abv in zip(recharge_list,rech_abv_list):
                        overall_count = overall_count + 1
                        rfd_top_com = recharge
                        name=str(overall_count) + '_var_' + str(var) + '_len_' + str(len_scale) + '_mean_' + str(mean) + '_seed_' + str(seed) + '_stor_' + str(storage) + '_rech_' + str(rech_abv)
                        output = "output"
                        dim_no = 2
                        parent_dir = CWD + '/setup'
                        # name of directory (entire path) of one single ogs run
                        dire = parent_dir + '/' + name
                        # make folders
                        if not os.path.exists(parent_dir):
                            os.mkdir(parent_dir)
                        if not os.path.exists(dire):
                            os.mkdir(dire)
                        pcs_type_flow = 'GROUNDWATER_FLOW'
                        var_name_flow = 'HEAD'
                        t_id = 'transect'
                        # ------------------------generate ogs base class----------------------------- #
                        ogs = OGS(task_root=dire+"/",
                                  task_id=t_id,
                                  output_dir=dire+"/"+output+"/")
                        # ------------------------  MSH -------------------------------- #
                        # generate a rectangular mesh in x-y-plane
                        ogs.msh.generate("rectangular", dim=2,
                                         mesh_origin=(0., 0.),
                                         element_no=(n_cellsx, n_cellsz),
                                         element_size=(1, 1))
                        # rotate mesh to obtain a cross section in x-z-plane
                        ogs.msh.rotate(angle=np.pi/2.0, rotation_axis=(1., 0., 0.))
                        # round nodes
                        ogs.msh.NODES[:, 1] = np.around(ogs.msh.NODES[:, 1], 0)
                        ogs.msh.NODES[:, 0] = np.around(ogs.msh.NODES[:, 0], 4)
                        ogs.msh.NODES[:, 2] = np.around(ogs.msh.NODES[:, 2], 4)

                        # ------------------------  heterogeneous field -------------------------------- #
                        # get the centroids of the mesh-elements to evaluate the field
                        cent = ogs.msh.centroids_flat
                        # split up in x and z coordinates
                        x = cent[:, 0]
                        z = cent[:, 2]
                        # generate the random field
                        cov_model = Gaussian(dim=dim, var=var, len_scale=len_scale)
                        srf = SRF(cov_model, mean=mean, seed=seed)
                        # use unstructured for a 2D vertical mesh
                        field = srf((x, z), mesh_type='unstructured', force_moments=True)#, mode_no=100)
                        # conductivities as log-normal distributed from the field data
                        cond = np.exp(field)
                        from scipy.stats.mstats import gmean, hmean
                        arimean = np.mean(cond)
                        harmean = hmean(cond)
                        geomean = gmean(cond)
                        print("The geometric mean is: " + str(geomean))
                        #plt.hist(field)
                        # show the heterogeneous field
                        plt.figure(figsize=(20, thickness/length*20))
                        cond_log = np.log10(cond)
                        plt.tricontourf(x, z, cond_log.T)
                        plt.colorbar(ticks=[np.min(cond_log),np.mean(cond_log),np.max(cond_log)])
                        plt.title("log-normal K field [log10 K]")
                        plt.savefig(dire + '/' + name + '.png', dpi=300, bbox_inches='tight')
                        plt.close()

                        # generate MPD file (media properties distributed)
                        mpd = MPD(ogs.task_id)
                        mpd.add_block(MSH_TYPE=pcs_type_flow,
                                      MMP_TYPE="PERMEABILITY",
                                      DIS_TYPE="ELEMENT",
                                      DATA=list(zip(range(len(cond)),cond))
                                      )

                        # add the field to the ogs model
                        ogs.add_mpd(mpd)
                        # export mesh and field for checking
                        ogs.msh.export_mesh(dire + '/' + t_id + "_hetero_field.vtk",
                                            file_format="vtk-ascii",
                                            add_data_by_id=cond)

                        # save a file with information about the generated field
                        field_info = open(dire+"/"+'field_info'+'.dat', 'w')
                        field_info.write('dim var len_scale mean seed geomean harmean arimean' + '\n' + (str(dim) + ' ' + str(var) + ' ' + str(len_scale) + ' ' + str(mean) + ' ' + str(seed)) + ' ' + str(geomean) + ' ' + str(harmean) + ' ' + str(arimean))
                        field_info.close()

                        # ------------------------  GLI -------------------------------- #
                        ogs.gli.add_points(points=[[0., 0., 0.],
                                                   [length, 0., 0.],
                                                   [length, 0., thickness],
                                                   [0., 0., thickness]],
                                           names=['A', 'B', 'C', 'D'])
                        # generate polylines from points as boundaries: always define polylines along positive x,y,z
                        ogs.gli.add_polyline(name='bottom', points=['A', 'B'])
                        ogs.gli.add_polyline(name='right', points=['B', 'C'])
                        ogs.gli.add_polyline(name='top', points=['D', 'C'])
                        ogs.gli.add_polyline(name='left', points=['A', 'D'])
                        # add the points and polylines based on the aquifer length
                        obs = []
                        percents_of_length = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.99, 1]
                        for percent in percents_of_length:
                            obs_loc = int(np.around(length*percent,0))
                            obs_str = 'obs_' + str(obs_loc).zfill(5)
                            obs.append(obs_str)
                            obs.sort()
                            ogs.gli.add_points(points=[obs_loc, 0., 0.], names=str('obs_' + str(obs_str) + '_bottom'))
                            ogs.gli.add_points(points=[obs_loc, 0., thickness], names=str('obs_' + str(obs_str) + '_top'))
                            ogs.gli.add_polyline(name=obs_str, points=[[obs_loc, 0., 0.], [obs_loc, 0., thickness]])

                        # --------------generate .rfd ------------------------- #
                        rfd_data = np.loadtxt(recharge)
                        # write array to .rfd, wenn mainkeyword als keyword in der funktion genutzt wird.
                        ogs.rfd.add_block(CURVES=rfd_data)
                        # kommentarzeile einfuegen
                        ogs.rfd.top_com = rfd_top_com
                        # --------------generate different ogs input classes------------------------- #

                        # --------------    BC  ------------------------- #
                        ogs.bc.add_block(PCS_TYPE=pcs_type_flow,
                                         PRIMARY_VARIABLE=var_name_flow,
                                         GEO_TYPE=[['POLYLINE', 'right']],
                                         DIS_TYPE=[['CONSTANT', thickness]])
                        # --------------    IC  ------------------------- #
                        ogs.ic.add_block(PCS_TYPE=pcs_type_flow,
                                         PRIMARY_VARIABLE=var_name_flow,
                                         GEO_TYPE='DOMAIN',
                                         DIS_TYPE=[['CONSTANT', thickness]])
                        # --------------    MFP ------------------------- #
                        ogs.mfp.add_block(FLUID_TYPE='WATER',
                                          DENSITY=[[1, 0.9997e+3]],
                                          VISCOSITY=[[1, 1.309e-3]])
                        # --------------    MMP ------------------------- #
                        ogs.mmp.add_block(GEOMETRY_DIMENSION=dim_no,
                                          STORAGE=[[1, storage]],
                                          PERMEABILITY_TENSOR=[['ISOTROPIC', 1]],
                                          PERMEABILITY_DISTRIBUTION=ogs.task_id+'.mpd',
                                          #POROSITY='0.35'
                                          )
                        # --------------    NUM ------------------------- #
                        ogs.num.add_block(PCS_TYPE=pcs_type_flow,
                                          # method error_tolerance max_iterations theta precond storage
                                          LINEAR_SOLVER=[[2, 1, 1.0e-10, 1000, 1.0, 100, 4]],
                                          ELE_GAUSS_POINTS=3,
                                          #NON_LINEAR_ITERATION=[['PICARD', 'ERNORM', 20, 0, 1e-6]]
                                          )
                        # --------------    OUT ------------------------- #
                        ogs.out.add_block(PCS_TYPE=pcs_type_flow,
                                          NOD_VALUES=[[var_name_flow],
                                                      ['VELOCITY_X1'],
                                                      ['VELOCITY_Z1']],
                                          ELE_VALUES=[['VELOCITY1_X'],
                                                      ['VELOCITY1_Z']],
                                          GEO_TYPE='DOMAIN',
                                          DAT_TYPE='PVD',
                                          TIM_TYPE=[['STEPS', 100]])

                        # set the output for every observation point
                        for obs_point in obs:
                            ogs.out.add_block(PCS_TYPE=pcs_type_flow,
                                              NOD_VALUES=[[var_name_flow],
                                                          ['VELOCITY_X1'],
                                                          ['VELOCITY_Z1']],
                                              GEO_TYPE=[['POLYLINE', obs_point]],
                                              DAT_TYPE='TECPLOT',
                                              TIM_TYPE=[['STEPS', 1]])

                        for state in ['steady','transient']:
                            # --------------    ST  ------------------------- #
                            if state == 'transient':
                                ogs.st.reset()
                                ogs.st.add_block(PCS_TYPE=pcs_type_flow,
                                                 PRIMARY_VARIABLE=var_name_flow,
                                                 GEO_TYPE=[['POLYLINE', 'top']],
                                                 DIS_TYPE=[['CONSTANT_NEUMANN', 1]],
                                                 TIM_TYPE=[['CURVE',1]])
                            if state == 'steady':
                                ogs.st.add_block(PCS_TYPE=pcs_type_flow,
                                                PRIMARY_VARIABLE=var_name_flow,
                                                GEO_TYPE=[['POLYLINE', 'top']],
                                                DIS_TYPE=[['CONSTANT_NEUMANN', np.mean(rfd_data[:,1])]])

                            # --------------    PCS ------------------------- #
                            if state == 'transient':
                                ogs.pcs.reset()
                                ogs.pcs.add_block(PCS_TYPE=pcs_type_flow,
                                                  NUM_TYPE='NEW',
                                                  PRIMARY_VARIABLE=var_name_flow,
                                                  RELOAD=[[2,1]],
                                                  BOUNDARY_CONDITION_OUTPUT=[[]])
                            if state == 'steady':
                                ogs.pcs.add_block(PCS_TYPE=pcs_type_flow,
                                                  TIM_TYPE='STEADY',
                                                  NUM_TYPE='NEW',
                                                  PRIMARY_VARIABLE=var_name_flow,
                                                  RELOAD=[[1,1]],
                                                  BOUNDARY_CONDITION_OUTPUT=[[]])
                            # --------------    TIM ------------------------- #
                            if state == 'transient':
                                ogs.tim.reset()
                                ogs.tim.add_block(PCS_TYPE=pcs_type_flow,
                                                  TIME_START=time_start,
                                                  TIME_END=time_end,
                                                  TIME_STEPS=zip(time_steps, step_size))

                            # --------------run OGS simulation------------------------------------------- #
                            ogs.write_input()
                            if state == 'steady':
                               file = open(dire+"/"+t_id+'.tim', 'w')
                               file.write('#STOP')
                               file.close()
                            if state == 'steady':
                                print("calculating steady state...")
                                ogs.run_model(ogs_root='/home/houben/OGS_source/ogs')



#print("run model")
#ogs.run_model(ogs_root='/Users/houben/PhD/ogs5/executable/ogs5_bugfix_RWPT')

#from ogs5py.reader import readtec_polyline
#tecs = readtec_polyline(task_id=name_of_project_ogs,task_root=path_to_project)

shutil.copyfile(str(CWD) + '/' + file_name, parent_dir + '/' + file_name)
print('OGS-model generation finished!')
