#! usr/bin/env python

import os, string

'''
In this script, you enter the parameters you want in the simulation, as well as 
the name of the executable, the job name, and the allocation.
The script then loops through the parameters you've given and generates the 
appropriate input file, then generates the submit script to run the code. 
Finally, it copies the executable to the appropriate directory.

It's up to you to run the scripts once they're generated.

'''

dim = 2 
dt = 0.05 
m = 1.
#list of values for chemical potential mu
mu_list = [-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0]
#list of values for rotational frequency wx
w_list = [0.0] 
#list of values for harmonic trap frequency wt
w_trap_list = [0.0] 
#list of values for interaction strength lambda
l_list = [0.0]
#list of values for number of spatial lattice sites
Nx_list = [21]
#list of values for number of temporal lattice site
Nt_list = [80]
#list of values for number of steps in the Langevin evolution
nL_list = [1000000]
#autocorrelation spacing
acf_spacing = 10
#thermalization step (-1 means we take samples without waiting to thermalize)
therm_step = -1
#list of epsilon values 
eps_list = [0.01]

script_name = "ComplexLangevin3d.gnu.haswell.TPROF.MPI.OMP.ex"
job_name = "amrex_cl_2D_free_gas"
allocation = "m3764" #this is our current allocation through the ERCAP grant
email = "caseyb@bu.edu"


def generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2,
	acf_spacing=1,num_plotfiles=1,therm_step=-1,
	use_hdf5=False,seed_init=8134,seed_run=61,max_grid_size=8):
	#calculate frequency of plotfile interval -- we want to have exactly num_plotfiles
	plot_int = int(nL/num_plotfiles)
	#generate file extension 
	file_ext = "d"+str(dim)+"_w"+str(w)+"_wt"+str(w_t)+"_lambda"+str(l)
	file_ext = file_ext+"_m"+str(m)+"_Nx"+str(Nx)+"_Nt"+str(Nt)+"_dt"+str(dt)
	file_ext = file_ext+"_nL"+str(nL)+"_eps"+str(eps)+"_mu"+str(mu)	
	#create directory
	currdir = os.getcwd()
	work_dir = currdir+'/nrrb_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	'''
		WRITE THE INPUT FILE
	'''
	filename = work_dir+"/inputs_"+file_ext
	input_file = open(filename,'w')
	input_file.write("# Physical Input Parameters\n")
	input_file.write("nrrb.m = "+str(m)+'\n')
	input_file.write("nrrb.l = "+str(l)+'\n')
	input_file.write("nrrb.w = "+str(w)+'\n')
	input_file.write("nrrb.w_t = "+str(w_t)+'\n')
	input_file.write("nrrb.dtau = "+str(dt)+'\n')
	input_file.write("nrrb.mu = "+str(mu)+'\n')
	input_file.write("nrrb.eps = "+str(eps)+'\n')
	input_file.write("nrrb.circulation_radius_1 = "+str(circ1)+'\n')
	input_file.write("nrrb.circulation_radius_2 = "+str(circ2)+'\n')
	input_file.write("\n# Set a seed for the random number generator\n")
	input_file.write("nrrb.seed_init = "+str(seed_init)+'\n')
	input_file.write("nrrb.seed_run = "+str(seed_run)+'\n')
	input_file.write("\n#Domain Size and Grid Decomposition\n")
	input_file.write("n_cell = "+str(Nx)+" "+str(Nx)+" "+str(Nt)+'\n')
	input_file.write("max_grid_size = "+str(max_grid_size)+'\n')
	input_file.write("\n# Set Fab tile size (default is 1024000 x 8 x 8 for 3D)\n")
	input_file.write("\n# Periodicity and Domain Boundary Conditions\n")
	'''
	For a given dimension, if is_periodic[i] = 1, then the lo/hi BC type 
	should be BCType::int_dir.
	
	Here are the supported boundary conditions in AMReX:
		BCType::reflect_odd  = -1
		BCType::int_dir      =  0
		BCType::reflect_even =  1
		BCType::foextrap     =  2
		BCType::ext_dir      =  3
		BCType::hoextrap     =  4
	'''
	input_file.write("is_periodic = 0 0 1 # Non-periodic in space, periodic in time\n")
	input_file.write("domain_lo_bc_types = 3 3 0 # External Dirichlet in space, interior BCs (periodic) in time\n")
	input_file.write("domain_hi_bc_types = 3 3 0 # External Dirichlet in space, interior BCs (periodic) in time\n")
	input_file.write("\n#Stopping criteria\n")
	input_file.write("nsteps = "+str(nL)+'\n')
	input_file.write("\n#Frequency for computing observables\n")
	input_file.write("autocorrelation_step = "+str(acf_spacing)+'\n')
	input_file.write("\n#I/O\n")
	input_file.write("# only enable plotting after taking a number of steps set by \"plot_after_step\"\n")
	input_file.write("plot_after_step = "+str(therm_step)+'\n')
	input_file.write("# plot every \"plot_int\" steps \n")
	input_file.write("plot_int = "+str(plot_int)+'\n')
	input_file.write("set logfile name \n")
	input_file.write("observable_log_file = \"observables_"+file_ext+".log\"\n")
	input_file.write("toggle HDF5 use \n")
	input_file.write("use_hdf5 = "+str(use_hdf5)+'\n')
	input_file.close()
	return file_ext

def copy_executable(script_name, file_ext):
	currdir = os.getcwd()
	work_dir = currdir+'/nrrb_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	source_path = currdir+"/"+script_name
	destination_path = work_dir+"/"
	copy_cmd = "cp "+source_path+" "+destination_path
	os.system(copy_cmd)

def generate_slurm_script(script_name,file_ext,job_name,allocation,num_nodes=2,
	tasks_per_node=2,cpus_per_task=32,queue="regular",walltime="48:00:00",
	email = "caseyb@bu.edu",omp_threads=16):
	#generates the sbatch file that you run with sbatch filename
	#should be paired with the appropriate input file somehow...
	filename = "cori.MPI.OMP.slurm"
	#create directory (if it doesn't already exist, which it should)
	currdir = os.getcwd()
	work_dir = currdir+'/nrrb_data_'+file_ext
	if not os.path.exists(work_dir):
		os.makedirs(work_dir)
	slurm_file = open(work_dir+'/'+filename,'w')
	slurm_file.write("#!/bin/bash\n#\n# This job requests Cori Haswell nodes\n\n")
	slurm_file.write("### For KNL, pass these options to --constraint: knl,quad,cache\n")
	slurm_file.write("#SBATCH --constraint=haswell\n#\n")
	slurm_file.write("# Number of nodes:\n#SBATCH --nodes="+str(num_nodes)+'\n')
	slurm_file.write("#\n# Assign 1 MPI task to each socket on the Haswell nodes:\n")
	slurm_file.write("#SBATCH --tasks-per-node="+str(tasks_per_node)+'\n')
	slurm_file.write("#\n#On Haswell, each socket has 32 CPUs (with hyperthreading) for 1 MPI task\n")
	slurm_file.write("#SBATCH --cpus-per-task="+str(cpus_per_task)+'\n')
	slurm_file.write("#\n# Which queue to run in: debug, regular, premium, etc. ...\n")
	slurm_file.write("#SBATCH --qos="+queue+'\n')
	slurm_file.write("#\n#Run for this much walltime: hh:mm:ss\n")
	slurm_file.write("#SBATCH --time="+walltime+'\n')
	slurm_file.write("#\n#Use this job name:\n")
	slurm_file.write("#SBATCH -J amrex_cl_test\n")
	slurm_file.write("#\n# Send notification emails here:\n")
	slurm_file.write("#SBATCH --mail-user="+email+"\n")
	slurm_file.write("#SBATCH --mail-type=ALL\n")
	slurm_file.write("#\n# Which allocation to use:\n")
	slurm_file.write("#SBATCH -A "+allocation+"\n")
	slurm_file.write("#\n# On the compute node, change to the directory we submitted from\n")
	slurm_file.write("cd $SLURM_SUBMIT_DIR\n\n")
	slurm_file.write("# OpenMP Configuration\n")
	slurm_file.write("export OMP_PLACES=threads\n")
	slurm_file.write("export OMP_PROC_BIND=true\n\n")
	slurm_file.write("## for Haswell nodes:\n")
	slurm_file.write("export OMP_NUM_THREADS="+str(omp_threads)+"\n")
	slurm_file.write("## for KNL nodes:\n")
	slurm_file.write("# export OMP_NUM_THREADS=68\n\n")
	slurm_file.write("srun --cpu_bind=cores ./"+script_name+" inputs_"+file_ext+"\n")
	slurm_file.close()



for Nx in Nx_list:
	circ1 = Nx/8
	circ2 = Nx/2 -1
	for Nt in Nt_list:
		for nL in nL_list:
			for eps in eps_list:
				for mu in mu_list:
					if not dim == 2:
						w = 0.0
						w_t = 0.0
						l = 0.0
						file_ext = generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2)
						generate_slurm_script(script_name,file_ext,job_name,allocation,email=email)
						copy_executable(script_name, file_ext)
					else:
						for l in l_list:
							for w_t in w_trap_list:
								for w in w_list:
									file_ext = generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2)
									generate_slurm_script(script_name,file_ext,job_name,allocation,email=email)
									copy_executable(script_name, file_ext)
