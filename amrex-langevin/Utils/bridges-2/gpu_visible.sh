#!/bin/bash
#
# This command assigns 1 GPU to each MPI rank.
#
# Uses a trick from https://medium.com/@jeffrey_91423/binding-to-the-right-gpu-in-mpi-cuda-programs-263ac753d232
#
# How to Use:
#
# - Interactive session on Bridges-2 GPU: get, e.g. 2 nodes like this
# -- `interact -p GPU -t 00:30:00 -N 2 --ntasks-per-node=8 -n 16 --gres=gpu:8`
# -- `mpirun ./gpu_visible.sh ./main.exe inputs` for an executable named "main.exe" that takes an "inputs" file.
#
# - Batch job on Bridges-2 GPU: put this into the batch script
# -- `mpirun ./gpu_visible.sh ./main.exe inputs` for an executable named "main.exe" that takes an "inputs" file.

# This sets the GPU ID that is visible to our MPI rank to be equal to its local node index
export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK

# This runs whatever commands follow the invocation of this script, e.g. to run the actual code
$@
