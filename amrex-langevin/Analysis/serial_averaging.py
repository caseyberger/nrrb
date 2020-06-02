import os
import numpy as np
from defs_cori import *


curr_dir = os.getcwd()+"/"
subdirectories = os.listdir(curr_dir)
data_directories = []
for subdir in subdirectories:
	if subdir.startswith("nrrb_data_"):
		work_dir = curr_dir+subdir
		if os.path.isdir(work_dir):
			all_files = os.listdir(work_dir)
			for filename in all_files:
				if filename.startswith("input"):
					input_filename = filename
				if filename.startswith("logfile"): 
					obs_filename = filename
			p = extract_parameters(work_dir,input_filename)
			Langevin_evol_data = get_raw_data(work_dir, obs_filename)
			Ntherm = int(0.2*len(Langevin_evol_data["Re n"]))
			average_observables(work_dir,Langevin_evol_data,p,Ntherm)
concatenate_obs_files(curr_dir)