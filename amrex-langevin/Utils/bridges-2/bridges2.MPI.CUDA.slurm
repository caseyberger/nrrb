#!/bin/bash

# SLURM script for PSC Bridges-2 GPU system.
# Requires the gpu_visible.sh script in the current directory.

# Sample SLURM command for Bridges-2 GPU:
#
# sbatch -p GPU -t 00:30:00 -N 1 --ntasks-per-node=8 --gres=gpu:8 bridges2.MPI.CUDA.slurm
#
# - Change walltime & number of nodes:
# -- `-t 00:30:00`: walltime in HH:MM:SS
# -- `-N 1`: number of nodes
#
# - Always use these options for the GPU nodes
# -- `-p GPU`: request GPU nodes
# -- `--ntasks-per-node=8`: 8 MPI ranks per node, 1 for each GPU
# -- `--gres=gpu:8`: use all 8 GPUs on each node

# On the compute node, change to the directory we submitted from
cd $SLURM_SUBMIT_DIR

# Run the executable using the gpu_visible wrapper to set 1 GPU per MPI rank
mpirun ./gpu_visible.sh ./main3d.gnu.MPI.CUDA.ex inputs
