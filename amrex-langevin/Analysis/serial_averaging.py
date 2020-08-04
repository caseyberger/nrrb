import os
import numpy as np
from defs_cori import *

curr_dir = os.getcwd()+"/"
subdirectories = os.listdir(curr_dir)
data_directories = []
for subdir in subdirectories:
    serial_averaging(subdir)
'''
	if subdir.startswith("nrrb_data_"):
		#print(subdir)
		work_dir = curr_dir+subdir
		if os.path.isdir(work_dir):
			all_files = os.listdir(work_dir)
			input_filename = "null"
			obs_filename = "null"
			for filename in all_files:
				if filename.startswith("input"):
					input_filename = filename
				if filename.startswith("logfile"): 
					obs_filename = filename
			if input_filename == "null":
				print("no input file")
				pass
			elif obs_filename == "null":
				print("no observables file")
				pass
			else:
				p = extract_parameters(work_dir,input_filename)
				Langevin_evol_data = get_raw_data(work_dir, obs_filename)
				Ntherm = int(0.2*float(p['nL']))
				average_observables(work_dir,Langevin_evol_data,p,Ntherm)
'''
concatenate_obs_files(curr_dir)
