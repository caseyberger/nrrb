#!/usr/bin/env python
import yt
import numpy as np
import os
import pandas as pd
from defs_cori import *

curr_dir = os.getcwd()+"/"
df = pd.DataFrame(columns=['step','nL','eps','Nx','Nt','dt','beta','mu','wz','wtr','lambda',
	'nRe','nIm','PhisqRe','PhisqIm','LzRe','LzIm'])
subdirectories = os.listdir(curr_dir)
for subdir in subdirectories:
    work_dir = curr_dir+subdir
    if os.path.isdir(work_dir):
        input_filename, plotfiles, steps, params = make_file_lists(subdir)
        nL_str = str(params["nL"])
        nL = params['nL']
        eps = params['eps']
        Nx = params['Nx']
        Nt = params['Nt']
        dt = params['dt']
        beta = params['beta']
        mu = params['mu']
        wz = params['wz']
        wtr = params['wtr']
        l = params['lambda']
        for step in np.arange(0,nL+1): #for step in steps:
            filename = "plt"+str(step).zfill(len(nL_str)-1)
            ds = yt.load(work_dir+"/"+filename)
            gd = GridData(ds)
            N = total_density(gd, Nx, Nt, beta, mu)
            Phisq = field_modulus_squared(gd, Nx, Nt)
            Lz = angular_momentum(gd, Nx, Nt)
            df = df.append({'step': step,  'nL': nL, 'eps': eps,'Nx': Nx, 'Nt': Nt,'dt': dt, 
            	'beta': beta, 'mu': mu, 'wz': wz, 'wtr': wtr, 'lambda': l,
            	'nRe': N.real, 'nIm': N.imag, 'PhisqRe': Phisq.real, 'PhisqIm': Phisq.imag,
            	'LzRe':Lz.real,'LzIm':Lz.imag }, ignore_index=True)
    else:
        pass
export_csv = df.to_csv(curr_dir+"raw_data.csv",index = None, header=True)