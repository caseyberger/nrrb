#! usr/bin/env python

import os, string
from analysis_defs import *

'''

This generates a bunch of input files and slurm scripts in directories 
named using the parameter values. All you need to do then is copy the 
main code file ("main3d.intel.haswell.MPI.OMP.ex") into each folder
and submit the batch script and it should run!

you can change the defaults in the input script and the slurm script if you want
'''

dim = 2
dt = 0.05
m = 1.
mu_list = [-0.2,-0.1,0.0]
w_list = [0.0]
w_trap_list = [0.0]
l_list = [0.0]
m_list = [1.0]
Nx_list = [41]
Nt_list = [160]
nL_list = [10000,100000,1000000]
tL_list = [100.]
script_name = "main3d.gnu.haswell.MPI.OMP.ex"
job_name = "amrex_cl_2D_tL_test"
allocation = "mp111"
acf_spacing = 10
use_hdf5 = True


for Nx in Nx_list:
	circ1 = Nx/8
	circ2 = Nx/2 -1
	for Nt in Nt_list:
		for nL in nL_list:
			nTherm = 0
			#nTherm = int(nL*0.2)
			for tL in tL_list:
				eps = float(tL)/float(nL)
				for mu in mu_list:
					if not dim == 2:
						w = 0.0
						w_t = 0.0
						l = 0.0
						file_ext = generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2)
						generate_slurm_script(script_name,file_ext,job_name,allocation)
						copy_executable(script_name, file_ext)
					else:
						for l in l_list:
							for w_t in w_trap_list:
								for w in w_list:
									file_ext = generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2,nTherm,acf_spacing=acf_spacing,num_plotfiles=nL/acf_spacing,use_hdf5=use_hdf5)
									generate_slurm_script(script_name,file_ext,job_name,allocation,use_hdf5=use_hdf5)
									#copy_executable(script_name, file_ext)
