import os
import numpy as np
import multiprocessing as mp
from defs_cori import *

def serial_averaging(subdir):
	empty_dirs = []
	curr_dir = os.getcwd()
	if subdir.startswith("nrrb_data_"):
		count = count + 1
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
				empty_dirs.append(subdir)
				pass
			elif obs_filename == "null":
				print("no observables file")
				empty_dirs.append(subdir)
				pass
			else:
				p = extract_parameters(work_dir,input_filename)
				Langevin_evol_data = get_raw_data(work_dir, obs_filename)
				Ntherm = int(0.2*len(Langevin_evol_data["Re n"]))
				average_observables(work_dir,Langevin_evol_data,p,Ntherm)
	return empty_dirs

#num_processors = mp.cpu_count()
num_processors = 4
pool = mp.Pool(num_processors)

curr_dir = os.getcwd()+"/"
subdirectories = os.listdir(curr_dir)
empty_dirs = pool.map(serial_averaging, [subdir for subdir in subdirectories])
concatenate_obs_files(curr_dir)
print(empty_dirs)
