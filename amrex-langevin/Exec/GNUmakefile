# CL_HOME defines the top level directory for the complex langevin code.
CL_HOME ?= ..

# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../Submodules/amrex

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_HDF5  ?= FALSE
COMP      = gnu
DIM       = 3

TINY_PROFILE  = TRUE
USE_PARTICLES = TRUE
PRECISION     = DOUBLE

Bpack   :=
Blocs   := .

USE_TEST_CONSTANT_RNG ?= FALSE
USE_TEST_SEED_RNG ?= FALSE
USE_TEST_UNIFORM_RNG ?= FALSE

include $(CL_HOME)/Make.CL
