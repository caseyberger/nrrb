#! usr/bin/env python

import os, string

def generate_input_file(dim,m,Nx,Nt,dt,nL,eps,mu,w_t,w,l,circ1,circ2,nTherm=0,acf_spacing=1,num_plotfiles=1,use_hdf5=True,max_grid_size=8,seed_init=8134,seed_run=61):
	#calculate frequency of plotfile interval -- we want to have exactly num_plotfiles
	plot_int = int(nL/num_plotfiles)
	#generate file extension 
	file_ext = "d"+str(dim)+"_w"+str(w)+"_wt"+str(w_t)+"_lambda"+str(l)+"_m"+str(m)+"_Nx"+str(Nx)+"_Nt"+str(Nt)+"_dt"+str(dt)+"_nL"+str(nL)+"_eps"+str(eps)+"_mu"+str(mu)	
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
	'''
	since this test inputs has max_grid_size = 8, using 4 x 4 tiles means
	e.g., an 8 x 8 x 8 grid will be tiled into 4 tiles of 8 x 4 x 4 so there 
	will be work for 4 OpenMP threads even if that 8 x 8 x 8 grid is the only
	grid owned by an MPI rank.
	'''
	input_file.write("fabarray.mfiter_tile_size = 1024000 4 4\n")
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
	input_file.write("plot_int = "+str(plot_int)+'\n')
	input_file.write("plot_after_step = "+str(nTherm)+'\n')
	if use_hdf5:
		input_file.write("use_hdf5=true\n")
	input_file.write("observable_log_file = \"observables_"+file_ext+".log\"\n")
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

def generate_slurm_script(script_name,file_ext,job_name,allocation,use_hdf5=True,num_nodes=2,tasks_per_node=2,cpus_per_task=32,queue="regular",walltime="48:00:00",email = "caseyb@bu.edu",omp_threads=16):
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
	slurm_file.write("# Load modules\n")
	slurm_file.write("module load gcc\n")
	if use_hdf5:
		slurm_file.write("module swap PrgEnv-intel PrgEnv-gnu\nmodule load cray-hdf5-parallel\n")
	slurm_file.write("\n")
	slurm_file.write("srun --cpu_bind=cores ./"+script_name+" inputs_"+file_ext+"\n")
	slurm_file.close()
