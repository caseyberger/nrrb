#!/bin/bash
#
# Number of nodes:
#SBATCH --nodes=32
#
#################
# KNL nodes
#################
#
# Requests Cori KNL nodes:
#SBATCH --constraint=knl,quad,cache
#
# Each KNL node has 68 physical cores
# Each physical core has 4 hardware threads
# This script will use 64 physical cores/KNL
# and 4 MPI tasks/KNL, so 16 physical cores/MPI task
#
# KNL: Assign 4 MPI tasks to each socket
#SBATCH --tasks-per-node=4
#
# If we're using 16 physical cores/MPI task,
# that's 16*4 = 64 hyperthreads = 64 "cpus"
#SBATCH --cpus-per-task=64
#
#################
# Queue & Job
#################
#
# Which queue to run in: debug, regular, premium, etc. ...
#SBATCH --qos=regular
#
# Run for this much walltime: hh:mm:ss
#SBATCH --time=00:30:00
#
# Use this job name:
#SBATCH -J cl_256_32
#
# Send notification emails here:
#SBATCH --mail-user=eugene.willcox@gmail.com
#SBATCH --mail-type=ALL
#
# Which allocation to use:
#SBATCH -A m3018

# On the compute node, change to the directory we submitted from
cd $SLURM_SUBMIT_DIR

# OpenMP Configuration
# This configuration ignores the hyperthreads
# and assigns 1 OpenMP thread/physical core
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=16

srun --cpu_bind=cores ./main3d.intel.mic-knl.MPI.OMP.ex inputs_scaling_256
