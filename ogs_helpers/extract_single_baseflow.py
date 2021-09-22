"""
Short script to extract the baseflow for an ogs model run.
"""

from calculate_flow import get_baseflow_from_polyline

task_id = "transect"
task_root = "/Users/houben/Desktop/eve_work/20191002_generate_ogs_homogeneous_baseflow_sa_2_left/1210_kf_2.927e-06_stor_0.0001_rech_whitenoise_FAILED"
single_file = task_root + "/" + "transect_ply_obs_01000_t8_GROUNDWATER_FLOW.tec"

get_baseflow_from_polyline(
    task_id, task_root, single_file, orientation="vertical", save_flow=True)
