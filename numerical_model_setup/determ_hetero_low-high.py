"""
This script generates OGS model setups for 2D, vertical and rectangular
model domain with mpi4py. Every model setup is first generated with a steady state
configuration, afterwards the OGS run starts, output is redirected into the
"steady" folder and steady input files are changed to input files with transient
settings.

How To
------

0) Install python 3.6 and all modules listed below.
    ogs5py
    sys
    os
    numpy
    matplotlib
    mpi4py
1) Set the ogs_root variable to the path/of/your/ogs/executable
2) Go through the code and change the parameters according to your need until
    the double ------ line
3) Run the script with the following command from your cmd-line.

    mpirun -n NUMBERofSLOTS python3 name_of_this_script.py the/working/directory NUMBERofSLOTS

    slots: number of cores which should be used

"""
# -------------------- import modules
import time
import sys
import os
import numpy as np
from ogs5py import OGS, MPD, MSH
from ogs5py.reader import readpvd, readtec_point
import matplotlib.pyplot as plt
from mpi4py import MPI

# -------------------- other configurations
ogs_root = "/home/houben/OGS_source/ogs"
# -------------------- get the arguments and pass them into variables
file_name = sys.argv[0]
CWD = sys.argv[1]
slots = int(sys.argv[2])

# -------------------- configure mpi4py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# -------------------- domain configuration (rectangular aquifer)
length = 1000
thickness = 30
# The number of cells in x and z direction (any integer).
n_cellsx = length
n_cellsz = thickness
# The size of cells in x and z direction (any float).
s_cellsx = 1
s_cellsz = 1

# -------------------- time configuration
# I prefer to work with seconds, because it is default.
time_start = 0
time_steps = np.array([365 * 30])
step_size = np.array([86400])
time_end = np.sum(time_steps * step_size)

# -------------------- parameter configuration
# Give a list of arbitrary length. Storage = specific storage
storage_list = [0.01, 0.0001]
# Give a list with strings to files which contain the recharge you want to apply.
recharge_path_list = [
    "/home/houben/recharge/recharge_daily.txt",
    "/home/houben/recharge/recharge_daily_30years_seconds_mm_mHM_estanis_danube.txt",
]
# According to your recharge_list give a name for each recharge.
rech_abv_list = ["whitenoise", "mHM"]
# Hydraulic conductivity values for each block. kf1 = left, kf2 = right
kf1 = 0.00001
kf2 = 0.001
# Set a start value for "overall_count" which is the index. I recommend to use
# as much digits as you will be generating new ogs models to end up with a
# consistant naming. I.e. if more than 100 take 1001 as start.
start = 1001
overall_count = start
# Generate a list where you want to set the boarders between high and low cond.
# zone.
border_list = np.arange(50, 900, 50).tolist() + [
    920,
    930,
    940,
    950,
    960,
    970,
    980,
    990,
]

# -------------------- model configurations
# Specify the PROCESS and PRIMARY_VARIABLE
pcs_type_flow = "GROUNDWATER_FLOW"
var_name_flow = "HEAD"
# Give it a name.
t_id = "transect"
# Specify the position of the observation points relative to the length of
# the aquifer. They will be rounded!!
percents_of_length = [
    0,
    0.05,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.85,
    0.9,
    0.92,
    0.94,
    0.96,
    0.99,
    1,
]

# -------------------- Make some folders
dim_no = 2
# Setups will be stored in CWD + "/setup", this folder is called
# the "parent directory"
parent_dir = CWD + "/setup"
# sleep a random amount of seconds so that folder creation has been completed
# for at least one core
time.sleep(np.random.rand()*5)
if not os.path.exists(parent_dir):
    os.mkdir(parent_dir)

# This loop iterates over all combinations of storage, recharge and border
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
for storage in storage_list:
    for recharge_path, rech_abv in zip(recharge_path_list, rech_abv_list):
        for border in border_list:
            # Only run the fllowing code "on the right rank":
            if (overall_count - start) % slots == rank:
                print(
                    "Rank "
                    + str(rank)
                    + " starts to generate the ogs setup files..."
                )
                # Use the path to the recharge.txt as top comment in the .rfd-file
                rfd_top_com = recharge_path
                # Name the folder
                name = (
                    str(overall_count)
                    + "_border_"
                    + str(border)
                    + "_stor_"
                    + str(storage)
                    + "_rech_"
                    + str(rech_abv)
                )
                # Name of directory (entire path) of one single ogs run
                dire = parent_dir + "/" + name
                # Make folders
                if not os.path.exists(dire):
                    os.mkdir(dire)
                # -------------------- generate ogs base class
                ogs = OGS(
                    task_root=dire + "/",
                    task_id=t_id,
                    output_dir=dire + "/" + "steady" + "/",
                )
                # -------------------- MSH
                # Generate a rectangular mesh in x-z-plane. The first and the
                # second part of the mesh will be generated independently.
                # Afterwards, the two parts will be combined to one single
                # mesh. First part on the left side. This generates a mesh in
                # x-y-plane which will be rotated in a few steps.
                ogs.msh.generate(
                    "rectangular",
                    dim=2,
                    mesh_origin=(0.0, 0.0),
                    element_no=(border, n_cellsz),
                    element_size=(s_cellsx, s_cellsz),
                )
                # Use a MATERIAL_ID to be able to distinguish afer combination.
                ogs.msh.MATERIAL_ID = 0
                # Second part on the right side.
                second_mesh = MSH()
                second_mesh.generate(
                    "rectangular",
                    dim=2,
                    mesh_origin=(border, 0.0),
                    element_no=(length - border, n_cellsz),
                    element_size=(s_cellsx, s_cellsz),
                )
                second_mesh.MATERIAL_ID = 1
                # Combine both meshes.
                ogs.msh.combine_mesh(second_mesh)
                # Rotate mesh to obtain a cross section in x-z-plane.
                ogs.msh.rotate(
                    angle=np.pi / 2.0, rotation_axis=(1.0, 0.0, 0.0)
                )
                # Round the values of the nodes.
                ogs.msh.NODES[:, 1] = np.around(ogs.msh.NODES[:, 1], 0)
                ogs.msh.NODES[:, 0] = np.around(ogs.msh.NODES[:, 0], 4)
                ogs.msh.NODES[:, 2] = np.around(ogs.msh.NODES[:, 2], 4)
                # Export the mesh if you like.
                # ogs.msh.export_mesh("Path/and/name/for/your/mesh.msh")

                # -------------------- GLI
                # Add points on every edge of the model and name them.
                ogs.gli.add_points(
                    points=[
                        [0.0, 0.0, 0.0],
                        [length, 0.0, 0.0],
                        [length, 0.0, thickness],
                        [0.0, 0.0, thickness],
                    ],
                    names=["A", "B", "C", "D"],
                )
                # Generate polylines from points as boundaries:
                # I always define polylines along positive x,y,z
                ogs.gli.add_polyline(name="bottom", points=["A", "B"])
                ogs.gli.add_polyline(name="right", points=["B", "C"])
                ogs.gli.add_polyline(name="top", points=["D", "C"])
                ogs.gli.add_polyline(name="left", points=["A", "D"])
                # Add the points and polylines based on the aquifer length and
                # desired relative position of the observation points.
                obs = []
                for percent in percents_of_length:
                    obs_loc = int(np.around(length * percent, 0))
                    # Name of the observation point.
                    obs_str = "obs_" + str(obs_loc).zfill(len(str(length)) + 1)
                    obs.append(obs_str)
                    obs.sort()
                    # Add upper and lower points:
                    ogs.gli.add_points(
                        points=[obs_loc, 0.0, 0.0],
                        names=str("obs_" + str(obs_str) + "_bottom"),
                    )
                    ogs.gli.add_points(
                        points=[obs_loc, 0.0, thickness],
                        names=str("obs_" + str(obs_str) + "_top"),
                    )
                    # Create a polyline in between.
                    ogs.gli.add_polyline(
                        name=obs_str,
                        points=[
                            [obs_loc, 0.0, 0.0],
                            [obs_loc, 0.0, thickness],
                        ],
                    )

                # -------------------- generate .rfd
                # Load the rfd data from the file.
                rfd_data = np.loadtxt(recharge_path)
                # Write array to .rfd.
                ogs.rfd.add_block(CURVES=rfd_data)
                # Insert a top comment.
                ogs.rfd.top_com = rfd_top_com

                # -------------------- ogs input classes
                # -------------------- BC
                # Set a constant head boundary on the right side.
                ogs.bc.add_block(
                    PCS_TYPE=pcs_type_flow,
                    PRIMARY_VARIABLE=var_name_flow,
                    GEO_TYPE=[["POLYLINE", "right"]],
                    DIS_TYPE=[["CONSTANT", thickness]],
                )

                # -------------------- IC
                # Set the initial condition to "fully saturated" = thickness.
                ogs.ic.add_block(
                    PCS_TYPE=pcs_type_flow,
                    PRIMARY_VARIABLE=var_name_flow,
                    GEO_TYPE="DOMAIN",
                    DIS_TYPE=[["CONSTANT", thickness]],
                )

                # -------------------- MFP
                # Standard fluid properties.
                ogs.mfp.add_block(
                    FLUID_TYPE="WATER",
                    DENSITY=[[1, 0.9997e3]],
                    VISCOSITY=[[1, 1.309e-3]],
                )

                # -------------------- MMP
                # First block for MATERIAL_ID = 0
                ogs.mmp.add_block(
                    GEOMETRY_DIMENSION=dim_no,
                    STORAGE=[[1, storage]],
                    PERMEABILITY_TENSOR=[["ISOTROPIC", kf1]],
                )
                # First block for MATERIAL_ID = 1
                ogs.mmp.add_block(
                    GEOMETRY_DIMENSION=dim_no,
                    STORAGE=[[1, storage]],
                    PERMEABILITY_TENSOR=[["ISOTROPIC", kf2]],
                )

                # -------------------- NUM
                # Set the linear solver.
                ogs.num.add_block(
                    PCS_TYPE=pcs_type_flow,
                    # method error_tolerance max_iterations theta precond storage
                    LINEAR_SOLVER=[[2, 1, 1.0e-10, 1000, 1.0, 100, 4]],
                    ELE_GAUSS_POINTS=3,
                )

                # -------------------- OUT
                # Write node and element velocities.
                ogs.out.add_block(
                    PCS_TYPE=pcs_type_flow,
                    NOD_VALUES=[
                        [var_name_flow],
                        ["VELOCITY_X1"],
                        ["VELOCITY_Y1"],
                        ["VELOCITY_Z1"],
                    ],
                    ELE_VALUES=[["VELOCITY1_X"],
                                ["VELOCITY1_Y"],
                                ["VELOCITY1_Z"]],
                    GEO_TYPE="DOMAIN",
                    # Write a PVD output for the whole domain.
                    DAT_TYPE="PVD",
                    # Write it for every 1000th time step.
                    TIM_TYPE=[["STEPS", 1000]],
                )

                # Get output for every observation point as tec file.
                for obs_point in obs:
                    ogs.out.add_block(
                        PCS_TYPE=pcs_type_flow,
                        NOD_VALUES=[
                            [var_name_flow],
                            ["VELOCITY_X1"],
                            ["VELOCITY_Y1"],
                            ["VELOCITY_Z1"],
                        ],
                        ELE_VALUES=[["VELOCITY1_X"],
                                    ["VELOCITY1_Y"],
                                    ["VELOCITY1_Z"]],
                        GEO_TYPE=[["POLYLINE", obs_point]],
                        DAT_TYPE="TECPLOT",
                        # For every time step.
                        TIM_TYPE=[["STEPS", 1]],
                    )
                # Save information about the block model in a file
                block_info = open(dire + "/" + "block_info" + ".dat", "w")
                block_info.write(
                    "border kf1 kf2"
                    + "\n"
                    + str(border)
                    + " "
                    + str(kf1)
                    + " "
                    + str(kf2)
                )
                block_info.close()

                # Write the steady state input files , run it , copy steady
                # files into "steady" folder and write the transient input
                # files.
                for state in ["steady", "transient"]:

                    # -------------------- ST
                    if state == "transient":
                        # Apply the recharge on the top of the model domain.
                        ogs.st.reset()
                        ogs.st.add_block(
                            PCS_TYPE=pcs_type_flow,
                            PRIMARY_VARIABLE=var_name_flow,
                            GEO_TYPE=[["POLYLINE", "top"]],
                            DIS_TYPE=[["CONSTANT_NEUMANN", 1]],
                            TIM_TYPE=[["CURVE", 1]],
                        )
                    if state == "steady":
                        ogs.st.add_block(
                            PCS_TYPE=pcs_type_flow,
                            PRIMARY_VARIABLE=var_name_flow,
                            GEO_TYPE=[["POLYLINE", "top"]],
                            DIS_TYPE=[
                                ["CONSTANT_NEUMANN", np.mean(rfd_data[:, 1])]
                            ],
                        )

                    # -------------------- PCS
                    if state == "transient":
                        ogs.pcs.reset()
                        ogs.pcs.add_block(
                            PCS_TYPE=pcs_type_flow,
                            NUM_TYPE="NEW",
                            PRIMARY_VARIABLE=var_name_flow,
                            RELOAD=[[2, 1]],
                            BOUNDARY_CONDITION_OUTPUT=[[]],
                        )
                    if state == "steady":
                        ogs.pcs.add_block(
                            PCS_TYPE=pcs_type_flow,
                            TIM_TYPE="STEADY",
                            NUM_TYPE="NEW",
                            PRIMARY_VARIABLE=var_name_flow,
                            RELOAD=[[1, 1]],
                            BOUNDARY_CONDITION_OUTPUT=[[]],
                        )

                    # -------------------- TIM
                    if state == "transient":
                        ogs.tim.reset()
                        ogs.tim.add_block(
                            PCS_TYPE=pcs_type_flow,
                            TIME_START=time_start,
                            TIME_END=time_end,
                            TIME_STEPS=zip(time_steps, step_size),
                        )

                    # -------------------- run OGS simulation
                    ogs.write_input()
                    if state == "steady":
                        file = open(dire + "/" + t_id + ".tim", "w")
                        file.write("#STOP")
                        file.close()
                        print("Running steady state for folder " + name + " on rank "  + str(rank))
                        ogs.run_model(ogs_root=ogs_root)
            # Increase the counter for the naming.
            # First folder will be equal to the value of start
            overall_count = overall_count + 1
