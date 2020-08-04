#!/bin/bash
#
# This job requests Cori Haswell nodes

### For KNL, pass these options to --constraint: knl,quad,cache
#SBATCH --constraint=haswell
#
# Number of nodes:
#SBATCH --nodes=2
#
# Assign 1 MPI task to each socket on the Haswell nodes:
#SBATCH --tasks-per-node=2
#
#On Haswell, each socket has 32 CPUs (with hyperthreading) for 1 MPI task
#SBATCH --cpus-per-task=32
#
# Which queue to run in: debug, regular, premium, etc. ...
#SBATCH --qos=regular
#
#Run for this much walltime: hh:mm:ss
#SBATCH --time=48:00:00
#
#Use this job name:
#SBATCH -J amrex_cl_test
#
# Send notification emails here:
#SBATCH --mail-user=caseyb@bu.edu
#SBATCH --mail-type=ALL
#
# Which allocation to use:
#SBATCH -A mp111
#
# On the compute node, change to the directory we submitted from
cd $SLURM_SUBMIT_DIR

# OpenMP Configuration
export OMP_PLACES=threads
export OMP_PROC_BIND=true

## for Haswell nodes:
export OMP_NUM_THREADS=16
## for KNL nodes:
# export OMP_NUM_THREADS=68

srun --cpu_bind=cores ./main3d.gnu.haswell.MPI.OMP.ex inputs_d2_w0.0_wt0.0_lambda0.0_m1.0_Nx21_Nt160_dt0.05_nL1000000_eps0.001_mu-0.1
