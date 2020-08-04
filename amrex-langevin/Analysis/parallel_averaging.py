import os
import numpy as np
import multiprocessing as mp
from defs_cori import *

'''
Last updated July 27, 2020
*pool.apply works better than pool.map, which was leaving out some directories (not entirely sure why)
'''

#num_processors = mp.cpu_count()
num_processors = 8
pool = mp.Pool(processes=num_processors)

curr_dir = os.getcwd()+"/"
subdirectories = os.listdir(curr_dir)
#count = pool.map(serial_averaging, [subdir for subdir in subdirectories])
count = [pool.apply(serial_averaging, args=(subdir,)) for subdir in subdirectories]
empty_dirs = concatenate_obs_files(curr_dir)