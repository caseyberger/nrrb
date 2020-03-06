#!/bin/bash
#
# Number of nodes:
#SBATCH --nodes=2
#
#################
# Haswell nodes
#################
#
# Requests Cori Haswell nodes:
#SBATCH --constraint=haswell
#
# Haswell: Assign 1 MPI task to each socket
#SBATCH --tasks-per-node=2
#
# Haswell: each socket has 32 CPUs (with hyperthreading)
#SBATCH --cpus-per-task=32
#
#################
# Queue & Job
#################
#
# Which queue to run in: debug, regular, premium, etc. ...
#SBATCH --qos=debug
#
# Run for this much walltime: hh:mm:ss
#SBATCH --time=00:30:00
#
# Use this job name:
#SBATCH -J amrex_cl
#
# Send notification emails here:
#SBATCH --mail-user=[EMAIL FOR NOTIFICATIONS]
#SBATCH --mail-type=ALL
#
# Which allocation to use:
#SBATCH -A m3018

# On the compute node, change to the directory we submitted from
cd $SLURM_SUBMIT_DIR

# OpenMP Configuration
export OMP_PLACES=threads
export OMP_PROC_BIND=true

## for Haswell nodes:
export OMP_NUM_THREADS=16

srun --cpu_bind=cores ./main3d.intel.haswell.MPI.OMP.ex inputs_mpi_omp_test
