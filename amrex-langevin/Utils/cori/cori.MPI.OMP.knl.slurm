#!/bin/bash
#
# Number of nodes:
#SBATCH --nodes=2
#
#################
# KNL nodes
#################
#
# Requests Cori KNL nodes:
#SBATCH --constraint=knl,quad,cache
#
# KNL: Assign 1 MPI task to each socket
#SBATCH --tasks-per-node=1
#
# KNL: each socket has 272 CPUs (with hyperthreading)
#SBATCH --cpus-per-task=272
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

# for KNL nodes:
export OMP_NUM_THREADS=68

srun --cpu_bind=cores ./main3d.intel.mic-knl.MPI.OMP.ex inputs_mpi_omp_test
