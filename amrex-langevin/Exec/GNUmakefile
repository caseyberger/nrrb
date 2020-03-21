# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex

DEBUG     = FALSE
USE_MPI   = FALSE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

USE_TEST_CONSTANT_RNG ?= FALSE
USE_TEST_SEED_RNG ?= FALSE
USE_TEST_UNIFORM_RNG ?= FALSE

ifeq ($(USE_TEST_CONSTANT_RNG), TRUE)
	DEFINES += -DTEST_CONSTANT_RNG
endif

ifeq ($(USE_TEST_SEED_RNG), TRUE)
	DEFINES += -DTEST_SEED_RNG
endif

ifeq ($(USE_TEST_UNIFORM_RNG), TRUE)
	DEFINES += -DTEST_UNIFORM_RNG
endif

include ../Source/Make.package

VPATH_LOCATIONS   += ../Source
INCLUDE_LOCATIONS += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
